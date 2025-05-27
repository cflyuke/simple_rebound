#include "simple_rebound.h"
#include <math.h>
#include <string.h>

// Gauss Radau spacings
static const double h[8] = { 
    0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 
    0.352624717113169637373907769648, 0.547153626330555383001448554766, 
    0.734210177215410531523210605558, 0.885320946839095768090359771030, 
    0.977520613561287501891174488626
};

// Other constants
static const double rr[28] = {
    0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 
    0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 
    0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 
    0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 
    0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 
    0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 
    0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 
    0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 
    0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 
    0.0921996667221917338008147
};

static const double c[21] = {
    -0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, 
    -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 
    0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, 
    -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, 
    -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 
    0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, 
    -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588
};

static const double d[21] = {
    0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 
    0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 
    0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 
    1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 
    0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 
    0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 
    0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588
};

// Safety factor for timestep control
static const double safety_factor = 0.25;

// Machine independent implementation of pow(*,1./7.)
static double sqrt7(double a) {
    double scale = 1;
    while(a < 1e-7 && isnormal(a)) {
        scale *= 0.1;
        a *= 1e7;
    }
    while(a > 1e2 && isnormal(a)) {
        scale *= 10;
        a *= 1e-7;
    }
    double x = 1.;
    for (int k = 0; k < 20; k++) {
        double x6 = x*x*x*x*x*x;
        x += (a/x6-x)/7.;
    }
    return x*scale;
}

static void free_dp7(IAS15_dp7* dp7) {
    free(dp7->p0); dp7->p0 = NULL;
    free(dp7->p1); dp7->p1 = NULL;
    free(dp7->p2); dp7->p2 = NULL;
    free(dp7->p3); dp7->p3 = NULL;
    free(dp7->p4); dp7->p4 = NULL;
    free(dp7->p5); dp7->p5 = NULL;
    free(dp7->p6); dp7->p6 = NULL;
}

static void clear_dp7(IAS15_dp7* dp7, int N3) {
    for (int k = 0; k < N3; k++) {
        dp7->p0[k] = 0.; dp7->p1[k] = 0.; dp7->p2[k] = 0.; dp7->p3[k] = 0.;
        dp7->p4[k] = 0.; dp7->p5[k] = 0.; dp7->p6[k] = 0.;
    }
}

static void realloc_dp7(IAS15_dp7* dp7, int N3) {
    dp7->p0 = realloc(dp7->p0, sizeof(double)*N3);
    dp7->p1 = realloc(dp7->p1, sizeof(double)*N3);
    dp7->p2 = realloc(dp7->p2, sizeof(double)*N3);
    dp7->p3 = realloc(dp7->p3, sizeof(double)*N3);
    dp7->p4 = realloc(dp7->p4, sizeof(double)*N3);
    dp7->p5 = realloc(dp7->p5, sizeof(double)*N3);
    dp7->p6 = realloc(dp7->p6, sizeof(double)*N3);
    clear_dp7(dp7, N3);
}

static inline void add_cs(double* p, double* csp, double inp) {
    const double y = inp - *csp;
    const double t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

static void copybuffers(IAS15_dp7* a, IAS15_dp7* b, int N3) {
    for (int i = 0; i < N3; i++) { 
        b->p0[i] = a->p0[i]; b->p1[i] = a->p1[i]; b->p2[i] = a->p2[i]; b->p3[i] = a->p3[i];
        b->p4[i] = a->p4[i]; b->p5[i] = a->p5[i]; b->p6[i] = a->p6[i];
    }
}

static void predict_next_step(double ratio, int N3, IAS15_dp7* _e, IAS15_dp7* _b, IAS15_dp7* e, IAS15_dp7* b) {
    if (ratio > 20.) {
        for(int k = 0; k < N3; ++k) {
            e->p0[k] = 0.; e->p1[k] = 0.; e->p2[k] = 0.; e->p3[k] = 0.; 
            e->p4[k] = 0.; e->p5[k] = 0.; e->p6[k] = 0.;
            b->p0[k] = 0.; b->p1[k] = 0.; b->p2[k] = 0.; b->p3[k] = 0.; 
            b->p4[k] = 0.; b->p5[k] = 0.; b->p6[k] = 0.;
        }
    } else {
        const double q1 = ratio;
        const double q2 = q1 * q1;
        const double q3 = q1 * q2;
        const double q4 = q2 * q2;
        const double q5 = q2 * q3;
        const double q6 = q3 * q3;
        const double q7 = q3 * q4;

        for(int k = 0; k < N3; ++k) {
            double be0 = _b->p0[k] - _e->p0[k];
            double be1 = _b->p1[k] - _e->p1[k];
            double be2 = _b->p2[k] - _e->p2[k];
            double be3 = _b->p3[k] - _e->p3[k];
            double be4 = _b->p4[k] - _e->p4[k];
            double be5 = _b->p5[k] - _e->p5[k];
            double be6 = _b->p6[k] - _e->p6[k];

            e->p0[k] = q1*(_b->p6[k]* 7.0 + _b->p5[k]* 6.0 + _b->p4[k]* 5.0 + _b->p3[k]* 4.0 + _b->p2[k]* 3.0 + _b->p1[k]*2.0 + _b->p0[k]);
            e->p1[k] = q2*(_b->p6[k]*21.0 + _b->p5[k]*15.0 + _b->p4[k]*10.0 + _b->p3[k]* 6.0 + _b->p2[k]* 3.0 + _b->p1[k]);
            e->p2[k] = q3*(_b->p6[k]*35.0 + _b->p5[k]*20.0 + _b->p4[k]*10.0 + _b->p3[k]* 4.0 + _b->p2[k]);
            e->p3[k] = q4*(_b->p6[k]*35.0 + _b->p5[k]*15.0 + _b->p4[k]* 5.0 + _b->p3[k]);
            e->p4[k] = q5*(_b->p6[k]*21.0 + _b->p5[k]* 6.0 + _b->p4[k]);
            e->p5[k] = q6*(_b->p6[k]* 7.0 + _b->p5[k]);
            e->p6[k] = q7* _b->p6[k];

            b->p0[k] = e->p0[k] + be0;
            b->p1[k] = e->p1[k] + be1;
            b->p2[k] = e->p2[k] + be2;
            b->p3[k] = e->p3[k] + be3;
            b->p4[k] = e->p4[k] + be4;
            b->p5[k] = e->p5[k] + be5;
            b->p6[k] = e->p6[k] + be6;
        }
    }
}

void ias15_alloc(Simulation* sim) {
    int N3 = 3 * sim->N;
    if (N3 > sim->ri_ias15.N_allocated) {
        realloc_dp7(&(sim->ri_ias15.g), N3);
        realloc_dp7(&(sim->ri_ias15.b), N3);
        realloc_dp7(&(sim->ri_ias15.csb), N3);
        realloc_dp7(&(sim->ri_ias15.e), N3);
        realloc_dp7(&(sim->ri_ias15.br), N3);
        realloc_dp7(&(sim->ri_ias15.er), N3);
        sim->ri_ias15.at = realloc(sim->ri_ias15.at, sizeof(double)*N3);
        sim->ri_ias15.x0 = realloc(sim->ri_ias15.x0, sizeof(double)*N3);
        sim->ri_ias15.v0 = realloc(sim->ri_ias15.v0, sizeof(double)*N3);
        sim->ri_ias15.a0 = realloc(sim->ri_ias15.a0, sizeof(double)*N3);
        sim->ri_ias15.csx = realloc(sim->ri_ias15.csx, sizeof(double)*N3);
        sim->ri_ias15.csv = realloc(sim->ri_ias15.csv, sizeof(double)*N3);
        sim->ri_ias15.csa0 = realloc(sim->ri_ias15.csa0, sizeof(double)*N3);
        
        for (int i = 0; i < N3; i++) {
            sim->ri_ias15.csx[i] = 0;
            sim->ri_ias15.csv[i] = 0;
        }
        sim->ri_ias15.N_allocated = N3;
    }
    if (sim->N > sim->ri_ias15.N_allocated_map) {
        sim->ri_ias15.map = realloc(sim->ri_ias15.map, sizeof(int)*sim->N);
        for (int i = 0; i < sim->N; i++) {
            sim->ri_ias15.map[i] = i;
        }
        sim->ri_ias15.N_allocated_map = sim->N;
    }
}

void ias15_reset(Simulation* sim) {
    sim->ri_ias15.N_allocated = 0;
    sim->ri_ias15.N_allocated_map = 0;
    free_dp7(&(sim->ri_ias15.g));
    free_dp7(&(sim->ri_ias15.e));
    free_dp7(&(sim->ri_ias15.b));
    free_dp7(&(sim->ri_ias15.csb));
    free_dp7(&(sim->ri_ias15.er));
    free_dp7(&(sim->ri_ias15.br));
    free(sim->ri_ias15.at); sim->ri_ias15.at = NULL;
    free(sim->ri_ias15.x0); sim->ri_ias15.x0 = NULL;
    free(sim->ri_ias15.v0); sim->ri_ias15.v0 = NULL;
    free(sim->ri_ias15.a0); sim->ri_ias15.a0 = NULL;
    free(sim->ri_ias15.csx); sim->ri_ias15.csx = NULL;
    free(sim->ri_ias15.csv); sim->ri_ias15.csv = NULL;
    free(sim->ri_ias15.csa0); sim->ri_ias15.csa0 = NULL;
    free(sim->ri_ias15.map); sim->ri_ias15.map = NULL;
}

int ias15_step(Simulation* sim) {
    ias15_alloc(sim);

    Particle* particles = sim->particles;
    int N = sim->N;
    int* map = sim->ri_ias15.map;
    const int N3 = 3*N;

    double* csx = sim->ri_ias15.csx; 
    double* csv = sim->ri_ias15.csv; 
    double* csa0 = sim->ri_ias15.csa0; 
    double* at = sim->ri_ias15.at; 
    double* x0 = sim->ri_ias15.x0; 
    double* v0 = sim->ri_ias15.v0; 
    double* a0 = sim->ri_ias15.a0; 
    IAS15_dp7* g = &(sim->ri_ias15.g);
    IAS15_dp7* e = &(sim->ri_ias15.e);
    IAS15_dp7* b = &(sim->ri_ias15.b);
    IAS15_dp7* csb = &(sim->ri_ias15.csb);
    IAS15_dp7* er = &(sim->ri_ias15.er);
    IAS15_dp7* br = &(sim->ri_ias15.br);
    
    for(int k = 0; k < N; k++) {
        int mk = map[k];
        x0[3*k]   = particles[mk].x;
        x0[3*k+1] = particles[mk].y;
        x0[3*k+2] = particles[mk].z;
        v0[3*k]   = particles[mk].vx;
        v0[3*k+1] = particles[mk].vy;
        v0[3*k+2] = particles[mk].vz;
        a0[3*k]   = particles[mk].ax;
        a0[3*k+1] = particles[mk].ay; 
        a0[3*k+2] = particles[mk].az;
    }
    
    for(int k = 0; k < N3; k++) {
        csa0[k] = 0;
    }
    
    for (int k = 0; k < N3; k++) {
        csb->p0[k] = 0.; csb->p1[k] = 0.; csb->p2[k] = 0.; csb->p3[k] = 0.;
        csb->p4[k] = 0.; csb->p5[k] = 0.; csb->p6[k] = 0.;
    }

    for(int k = 0; k < N3; k++) {
        g->p0[k] = b->p6[k]*d[15] + b->p5[k]*d[10] + b->p4[k]*d[6] + b->p3[k]*d[3]  + b->p2[k]*d[1]  + b->p1[k]*d[0]  + b->p0[k];
        g->p1[k] = b->p6[k]*d[16] + b->p5[k]*d[11] + b->p4[k]*d[7] + b->p3[k]*d[4]  + b->p2[k]*d[2]  + b->p1[k];
        g->p2[k] = b->p6[k]*d[17] + b->p5[k]*d[12] + b->p4[k]*d[8] + b->p3[k]*d[5]  + b->p2[k];
        g->p3[k] = b->p6[k]*d[18] + b->p5[k]*d[13] + b->p4[k]*d[9] + b->p3[k];
        g->p4[k] = b->p6[k]*d[19] + b->p5[k]*d[14] + b->p4[k];
        g->p5[k] = b->p6[k]*d[20] + b->p5[k];
        g->p6[k] = b->p6[k];
    }

    double t_beginning = sim->t;
    double predictor_corrector_error = 1e300;
    double predictor_corrector_error_last = 2;
    int iterations = 0; 
    
    while(1) {
        if(predictor_corrector_error < 1e-16) {
            break;
        }
        if(iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error) {
            break;
        }
        if (iterations >= 12) {
            sim->ri_ias15.iterations_max_exceeded++;
            break;
        }
        predictor_corrector_error_last = predictor_corrector_error;
        predictor_corrector_error = 0;
        iterations++;

        for(int n = 1; n < 8; n++) {
            sim->t = t_beginning + sim->dt * h[n];

            for(int i = 0; i < N; i++) {
                int mi = map[i];
                const int k0 = 3*i+0;
                const int k1 = 3*i+1;
                const int k2 = 3*i+2;

                double xk0, xk1, xk2;
                xk0 = -csx[k0] + ((((((((b->p6[k0]*7.*h[n]/9. + b->p5[k0])*3.*h[n]/4. + b->p4[k0])*5.*h[n]/7. + b->p3[k0])*2.*h[n]/3. + b->p2[k0])*3.*h[n]/5. + b->p1[k0])*h[n]/2. + b->p0[k0])*h[n]/3. + a0[k0])*sim->dt*h[n]/2. + v0[k0])*sim->dt*h[n];
                xk1 = -csx[k1] + ((((((((b->p6[k1]*7.*h[n]/9. + b->p5[k1])*3.*h[n]/4. + b->p4[k1])*5.*h[n]/7. + b->p3[k1])*2.*h[n]/3. + b->p2[k1])*3.*h[n]/5. + b->p1[k1])*h[n]/2. + b->p0[k1])*h[n]/3. + a0[k1])*sim->dt*h[n]/2. + v0[k1])*sim->dt*h[n];
                xk2 = -csx[k2] + ((((((((b->p6[k2]*7.*h[n]/9. + b->p5[k2])*3.*h[n]/4. + b->p4[k2])*5.*h[n]/7. + b->p3[k2])*2.*h[n]/3. + b->p2[k2])*3.*h[n]/5. + b->p1[k2])*h[n]/2. + b->p0[k2])*h[n]/3. + a0[k2])*sim->dt*h[n]/2. + v0[k2])*sim->dt*h[n];
                particles[mi].x = xk0 + x0[k0];
                particles[mi].y = xk1 + x0[k1];
                particles[mi].z = xk2 + x0[k2];
            }

            calculate_gravity(sim);

            for(int k = 0; k < N; ++k) {
                int mk = map[k];
                at[3*k]   = particles[mk].ax;
                at[3*k+1] = particles[mk].ay;  
                at[3*k+2] = particles[mk].az;
            }
            
            switch (n) {
                case 1: 
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p0[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p0[k] = gk/rr[0];
                        add_cs(&(b->p0[k]), &(csb->p0[k]), g->p0[k]-tmp);
                    } break;
                case 2: 
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p1[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p1[k] = (gk/rr[1] - g->p0[k])/rr[2];
                        tmp = g->p1[k] - tmp;
                        add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[0]);
                        add_cs(&(b->p1[k]), &(csb->p1[k]), tmp);
                    } break;
                case 3: 
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p2[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p2[k] = ((gk/rr[3] - g->p0[k])/rr[4] - g->p1[k])/rr[5];
                        tmp = g->p2[k] - tmp;
                        add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[1]);
                        add_cs(&(b->p1[k]), &(csb->p1[k]), tmp * c[2]);
                        add_cs(&(b->p2[k]), &(csb->p2[k]), tmp);
                    } break;
                case 4:
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p3[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p3[k] = (((gk/rr[6] - g->p0[k])/rr[7] - g->p1[k])/rr[8] - g->p2[k])/rr[9];
                        tmp = g->p3[k] - tmp;
                        add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[3]);
                        add_cs(&(b->p1[k]), &(csb->p1[k]), tmp * c[4]);
                        add_cs(&(b->p2[k]), &(csb->p2[k]), tmp * c[5]);
                        add_cs(&(b->p3[k]), &(csb->p3[k]), tmp);
                    } break;
                case 5:
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p4[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p4[k] = ((((gk/rr[10] - g->p0[k])/rr[11] - g->p1[k])/rr[12] - g->p2[k])/rr[13] - g->p3[k])/rr[14];
                        tmp = g->p4[k] - tmp;
                        add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[6]);
                        add_cs(&(b->p1[k]), &(csb->p1[k]), tmp * c[7]);
                        add_cs(&(b->p2[k]), &(csb->p2[k]), tmp * c[8]);
                        add_cs(&(b->p3[k]), &(csb->p3[k]), tmp * c[9]);
                        add_cs(&(b->p4[k]), &(csb->p4[k]), tmp);
                    } break;
                case 6:
                    for(int k = 0; k < N3; ++k) {
                        double tmp = g->p5[k];
                        double gk = at[k];
                        double gk_cs = 0.0;
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g->p5[k] = (((((gk/rr[15] - g->p0[k])/rr[16] - g->p1[k])/rr[17] - g->p2[k])/rr[18] - g->p3[k])/rr[19] - g->p4[k])/rr[20];
                        tmp = g->p5[k] - tmp;
                        add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[10]);
                        add_cs(&(b->p1[k]), &(csb->p1[k]), tmp * c[11]);
                        add_cs(&(b->p2[k]), &(csb->p2[k]), tmp * c[12]);
                        add_cs(&(b->p3[k]), &(csb->p3[k]), tmp * c[13]);
                        add_cs(&(b->p4[k]), &(csb->p4[k]), tmp * c[14]);
                        add_cs(&(b->p5[k]), &(csb->p5[k]), tmp);
                    } break;
                case 7:
                    {
                        double maxak = 0.0;
                        double maxb6ktmp = 0.0;
                        for(int k = 0; k < N3; ++k) {
                            double tmp = g->p6[k];
                            double gk = at[k];
                            double gk_cs = 0.0;
                            add_cs(&gk, &gk_cs, -a0[k]);
                            add_cs(&gk, &gk_cs, csa0[k]);
                            g->p6[k] = ((((((gk/rr[21] - g->p0[k])/rr[22] - g->p1[k])/rr[23] - g->p2[k])/rr[24] - g->p3[k])/rr[25] - g->p4[k])/rr[26] - g->p5[k])/rr[27];
                            tmp = g->p6[k] - tmp;    
                            add_cs(&(b->p0[k]), &(csb->p0[k]), tmp * c[15]);
                            add_cs(&(b->p1[k]), &(csb->p1[k]), tmp * c[16]);
                            add_cs(&(b->p2[k]), &(csb->p2[k]), tmp * c[17]);
                            add_cs(&(b->p3[k]), &(csb->p3[k]), tmp * c[18]);
                            add_cs(&(b->p4[k]), &(csb->p4[k]), tmp * c[19]);
                            add_cs(&(b->p5[k]), &(csb->p5[k]), tmp * c[20]);
                            add_cs(&(b->p6[k]), &(csb->p6[k]), tmp);

                            // 监控b.p6[k]相对于at[k]的变化，用于自适应步长控制
                            if (sim->ri_ias15.adaptive_mode != 0) {
                                const double ak = fabs(at[k]);
                                if (isnormal(ak) && ak > maxak) {
                                    maxak = ak;
                                }
                                const double b6ktmp = fabs(tmp);
                                if (isnormal(b6ktmp) && b6ktmp > maxb6ktmp) {
                                    maxb6ktmp = b6ktmp;
                                }
                            } else {
                                // 对于adaptive_mode == 0，使用传统的逐粒子误差估计作为收敛判据
                                const double ak = fabs(at[k]);
                                const double b6ktmp = fabs(tmp);
                                if (isnormal(ak) && ak > 0) {
                                    const double errork = b6ktmp/ak;
                                    if (isnormal(errork) && errork > predictor_corrector_error) {
                                        predictor_corrector_error = errork;
                                    }
                                }
                            }
                        }
                        if (sim->ri_ias15.adaptive_mode != 0) {
                            if (maxak > 0 && isnormal(maxak)) {
                                predictor_corrector_error = maxb6ktmp/maxak;
                            }
                        }
                        break;
                    }
            }
        }
    }
    
    // Set time back to initial value
    sim->t = t_beginning;
    
    // Calculate new timestep using adaptive control
    const double dt_done = sim->dt;
    double dt_new = dt_done;
    
    // Skip adaptive timestep control if exact finish time is required
    if (sim->ri_ias15.epsilon > 0 && !sim->exact_finish_time) {
        if (sim->ri_ias15.adaptive_mode == 2) {
            // PRS23 adaptive timestep control (Pham, Rein, Spiegel 2023)
            double min_timescale2 = INFINITY;
            unsigned int Nreal = N; // 实际粒子数量，不包括变分粒子
            
            for(unsigned int i = 0; i < Nreal; i++) {
                double a0i = 0; // 初始加速度平方
                double y2 = 0;  // 末端加速度平方
                double y3 = 0;  // (jerk * dt_done)^2
                double y4 = 0;  // (snap * dt_done^2)^2
                
                for(unsigned int k = 3*i; k < 3*(i+1); k++) {
                    a0i += a0[k]*a0[k];
                    // 计算末端加速度
                    double tmp = a0[k] + b->p0[k] + b->p1[k] + b->p2[k] + b->p3[k] + b->p4[k] + b->p5[k] + b->p6[k];
                    y2 += tmp*tmp;
                    // 计算jerk * dt_done
                    tmp = b->p0[k] + 2.* b->p1[k] + 3.* b->p2[k] + 4.* b->p3[k] + 5.* b->p4[k] + 6.* b->p5[k] + 7.* b->p6[k];
                    y3 += tmp*tmp;
                    // 计算snap * dt_done^2
                    tmp = 2.* b->p1[k] + 6.* b->p2[k] + 12.* b->p3[k] + 20.* b->p4[k] + 30.* b->p5[k] + 42.* b->p6[k];
                    y4 += tmp*tmp;
                }
                
                if (!isnormal(a0i)) {
                    // 跳过没有加速度或加速度为无穷大/NaN的粒子
                    continue;
                }
                
                double timescale2 = 2.*y2/(y3+sqrt(y4*y2)); // PRS23公式
                if (isnormal(timescale2) && timescale2 < min_timescale2) {
                    min_timescale2 = timescale2;
                }
            }
            
            if (isnormal(min_timescale2)) {
                // 数值因子用于匹配默认epsilon的时间步长
                dt_new = sqrt(min_timescale2) * dt_done * sqrt7(sim->ri_ias15.epsilon*5040.0);
            } else {
                dt_new = dt_done / safety_factor; // 默认情况下稍微增加时间步长
            }
        } else {
            // 传统的自适应时间步长控制
            double integrator_error = predictor_corrector_error;
            
            if (isnormal(integrator_error)) {
                dt_new = sqrt7(sim->ri_ias15.epsilon/integrator_error) * dt_done;
            } else {
                dt_new = dt_done / safety_factor; // 增加时间步长
            }
        }
        
        // 应用最小时间步长约束
        if (fabs(dt_new) < sim->ri_ias15.min_dt) {
            dt_new = copysign(sim->ri_ias15.min_dt, dt_new);
        }
        
        // 检查是否需要显著减少时间步长
        if (fabs(dt_new/dt_done) < safety_factor) {
            // 重置粒子到初始状态并用更小的时间步长重试
            for(int k = 0; k < N; ++k) {
                int mk = map[k];
                particles[mk].x = x0[3*k+0];
                particles[mk].y = x0[3*k+1];
                particles[mk].z = x0[3*k+2];
                particles[mk].vx = v0[3*k+0];
                particles[mk].vy = v0[3*k+1];
                particles[mk].vz = v0[3*k+2];
                particles[mk].ax = a0[3*k+0];
                particles[mk].ay = a0[3*k+1];
                particles[mk].az = a0[3*k+2];
            }
            sim->dt = dt_new;
            
            // 如果不是第一步，预测下一步系数
            if (sim->ri_ias15.dt_last_done != 0.0) {
                double ratio = sim->dt / sim->ri_ias15.dt_last_done;
                predict_next_step(ratio, N3, er, br, e, b);
            }
            
            return 0; // 步长被拒绝，重试
        }
        
        // 限制时间步长增加
        if (fabs(dt_new/dt_done) > 1.0) {
            if (dt_new/dt_done > 1.0/safety_factor) {
                dt_new = dt_done / safety_factor;
            }
        }
        
        sim->dt = dt_new;
    }

    // Find new position and velocity values at end of the sequence
    for(int k = 0; k < N3; ++k) {
        add_cs(&(x0[k]), &(csx[k]), b->p6[k]/72.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p5[k]/56.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p4[k]/42.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p3[k]/30.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p2[k]/20.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p1[k]/12.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), b->p0[k]/6.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), a0[k]/2.*dt_done*dt_done);
        add_cs(&(x0[k]), &(csx[k]), v0[k]*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p6[k]/8.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p5[k]/7.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p4[k]/6.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p3[k]/5.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p2[k]/4.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p1[k]/3.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), b->p0[k]/2.*dt_done);
        add_cs(&(v0[k]), &(csv[k]), a0[k]*dt_done);
    }

    sim->t += dt_done;
    sim->ri_ias15.dt_last_done = dt_done;

    // Update particle positions and velocities
    for(int k = 0; k < N; ++k) {
        int mk = map[k];
        particles[mk].x = x0[3*k+0];
        particles[mk].y = x0[3*k+1];
        particles[mk].z = x0[3*k+2];
        particles[mk].vx = v0[3*k+0];
        particles[mk].vy = v0[3*k+1];
        particles[mk].vz = v0[3*k+2];
    }
    
    // Copy buffers and predict next step
    copybuffers(e, er, N3);       
    copybuffers(b, br, N3);       
    double ratio = sim->dt/dt_done;
    predict_next_step(ratio, N3, e, b, e, b);
    
    return 1; // Success
}

// Main IAS15 integrator function with adaptive timestep
void integrator_ias15(Simulation* sim) {
    if (!sim) return;
    
    // Initialize IAS15 parameters if not set
    if (sim->ri_ias15.epsilon == 0.0) {
        sim->ri_ias15.epsilon = 1e-9;  // Default tolerance
    }
    if (sim->ri_ias15.min_dt == 0.0) {
        sim->ri_ias15.min_dt = 1e-12;  // Default minimum timestep
    }
    sim->ri_ias15.adaptive_mode = 2;  // Use PRS23 adaptive mode
    
    // Calculate initial accelerations
    calculate_gravity(sim);
    
    // Try until a step was successful
    while(!ias15_step(sim));
}
