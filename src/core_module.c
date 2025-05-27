#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "simple_rebound.h"

// Python对象包装C结构
typedef struct {
    PyObject_HEAD
    Simulation* sim;
} SimulationObject;

// 前向声明
static PyTypeObject SimulationType;

// 创建新的模拟对象
static PyObject* Simulation_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    SimulationObject *self;
    self = (SimulationObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->sim = sim_create();
        if (self->sim == NULL) {
            Py_DECREF(self);
            PyErr_SetString(PyExc_RuntimeError, "Failed to create simulation");
            return NULL;
        }
    }
    return (PyObject *)self;
}

// 销毁模拟对象
static void Simulation_dealloc(SimulationObject *self) {
    if (self->sim) {
        sim_destroy(self->sim);
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

// 添加粒子方法
static PyObject* Simulation_add_particle(SimulationObject *self, PyObject *args, PyObject *kwargs) {
    double m = 0.0, x = 0.0, y = 0.0, z = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0, r = 0.0;
    
    static char *kwlist[] = {"m", "x", "y", "z", "vx", "vy", "vz", "r", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|dddddddd", kwlist,
                                     &m, &x, &y, &z, &vx, &vy, &vz, &r)) {
        return NULL;
    }
    
    sim_add_particle(self->sim, m, x, y, z, vx, vy, vz, r);
    Py_RETURN_NONE;
}

// 积分一步
static PyObject* Simulation_step(SimulationObject *self, PyObject *args) {
    integrator_step(self->sim);
    if (self->sim->collision != COLLISION_NONE) {
        check_collisions(self->sim);
    }
    Py_RETURN_NONE;
}

// 积分到指定时间
static PyObject* Simulation_integrate(SimulationObject *self, PyObject *args) {
    double tmax;
    if (!PyArg_ParseTuple(args, "d", &tmax)) {
        return NULL;
    }
    
    while (self->sim->t < tmax) {
        integrator_step(self->sim);
        if (self->sim->collision != COLLISION_NONE) {
            check_collisions(self->sim);
        }
    }
    
    Py_RETURN_NONE;
}

// 获取粒子数量
static PyObject* Simulation_get_N(SimulationObject *self, void *closure) {
    return PyLong_FromLong(self->sim->N);
}

// 获取时间
static PyObject* Simulation_get_t(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->t);
}

// 设置时间
static int Simulation_set_t(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete time attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "Time must be a float");
        return -1;
    }
    
    self->sim->t = PyFloat_AsDouble(value);
    return 0;
}

// 获取时间步长
static PyObject* Simulation_get_dt(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->dt);
}

// 设置时间步长
static int Simulation_set_dt(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete dt attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "dt must be a float");
        return -1;
    }
    
    self->sim->dt = PyFloat_AsDouble(value);
    return 0;
}

// 获取引力常数
static PyObject* Simulation_get_G(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->G);
}

// 设置引力常数
static int Simulation_set_G(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete G attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "G must be a float");
        return -1;
    }
    
    self->sim->G = PyFloat_AsDouble(value);
    return 0;
}

// 获取软化系数
static PyObject* Simulation_get_softening(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->softening);
}

// 设置软化系数
static int Simulation_set_softening(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete softening attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "softening must be a float");
        return -1;
    }
    
    self->sim->softening = PyFloat_AsDouble(value);
    return 0;
}

// 获取theta参数
static PyObject* Simulation_get_theta(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->theta);
}

// 设置theta参数
static int Simulation_set_theta(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete theta attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "theta must be a float");
        return -1;
    }
    
    self->sim->theta = PyFloat_AsDouble(value);
    return 0;
}

// 获取积分器类型
static PyObject* Simulation_get_integrator(SimulationObject *self, void *closure) {
    switch (self->sim->integrator) {
        case INTEGRATOR_LEAPFROG:
            return PyUnicode_FromString("leapfrog");
        case INTEGRATOR_IAS15:
            return PyUnicode_FromString("ias15");
        case INTEGRATOR_WHFAST:
            return PyUnicode_FromString("whfast");
        case INTEGRATOR_MERCURIUS:
            return PyUnicode_FromString("mercurius");
        default:
            return PyUnicode_FromString("unknown");
    }
}

// 设置积分器类型
static int Simulation_set_integrator(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete integrator attribute");
        return -1;
    }
    
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "integrator must be a string");
        return -1;
    }
    
    const char *integrator_str = PyUnicode_AsUTF8(value);
    if (strcmp(integrator_str, "leapfrog") == 0) {
        self->sim->integrator = INTEGRATOR_LEAPFROG;
    } else if (strcmp(integrator_str, "ias15") == 0) {
        self->sim->integrator = INTEGRATOR_IAS15;
    } else if (strcmp(integrator_str, "whfast") == 0) {
        self->sim->integrator = INTEGRATOR_WHFAST;
    } else if (strcmp(integrator_str, "mercurius") == 0) {
        self->sim->integrator = INTEGRATOR_MERCURIUS;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown integrator type");
        return -1;
    }
    
    return 0;
}

// 获取碰撞类型
static PyObject* Simulation_get_collision(SimulationObject *self, void *closure) {
    switch (self->sim->collision) {
        case COLLISION_NONE:
            return PyUnicode_FromString("none");
        case COLLISION_MERGE:
            return PyUnicode_FromString("merge");
        case COLLISION_BOUNCE:
            return PyUnicode_FromString("bounce");
        default:
            return PyUnicode_FromString("unknown");
    }
}

// 设置碰撞类型
static int Simulation_set_collision(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete collision attribute");
        return -1;
    }
    
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "collision must be a string");
        return -1;
    }
    
    const char *collision_str = PyUnicode_AsUTF8(value);
    if (strcmp(collision_str, "none") == 0) {
        self->sim->collision = COLLISION_NONE;
    } else if (strcmp(collision_str, "merge") == 0) {
        self->sim->collision = COLLISION_MERGE;
    } else if (strcmp(collision_str, "bounce") == 0) {
        self->sim->collision = COLLISION_BOUNCE;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown collision type");
        return -1;
    }
    
    return 0;
}

// 获取引力计算方法
static PyObject* Simulation_get_gravity_method(SimulationObject *self, void *closure) {
    switch (self->sim->gravity_method) {
        case GRAVITY_DIRECT:
            return PyUnicode_FromString("direct");
        case GRAVITY_TREE:
            return PyUnicode_FromString("tree");
        default:
            return PyUnicode_FromString("unknown");
    }
}

// 设置引力计算方法
static int Simulation_set_gravity_method(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete gravity_method attribute");
        return -1;
    }
    
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "gravity_method must be a string");
        return -1;
    }
    
    const char *method_str = PyUnicode_AsUTF8(value);
    if (strcmp(method_str, "direct") == 0) {
        self->sim->gravity_method = GRAVITY_DIRECT;
    } else if (strcmp(method_str, "tree") == 0) {
        self->sim->gravity_method = GRAVITY_TREE;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown gravity method type");
        return -1;
    }
    
    return 0;
}

// 获取碰撞检测方法
static PyObject* Simulation_get_collision_detection(SimulationObject *self, void *closure) {
    switch (self->sim->collision_detection) {
        case COLLISION_DETECTION_DIRECT:
            return PyUnicode_FromString("direct");
        case COLLISION_DETECTION_SPATIAL:
            return PyUnicode_FromString("spatial");
        default:
            return PyUnicode_FromString("unknown");
    }
}

// 设置碰撞检测方法
static int Simulation_set_collision_detection(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete collision_detection attribute");
        return -1;
    }
    
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "collision_detection must be a string");
        return -1;
    }
    
    const char *method_str = PyUnicode_AsUTF8(value);
    if (strcmp(method_str, "direct") == 0) {
        self->sim->collision_detection = COLLISION_DETECTION_DIRECT;
    } else if (strcmp(method_str, "spatial") == 0) {
        self->sim->collision_detection = COLLISION_DETECTION_SPATIAL;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown collision detection method type");
        return -1;
    }
    
    return 0;
}

// 获取边界类型
static PyObject* Simulation_get_boundary(SimulationObject *self, void *closure) {
    switch (self->sim->boundary) {
        case BOUNDARY_OPEN:
            return PyUnicode_FromString("open");
        case BOUNDARY_PERIODIC:
            return PyUnicode_FromString("periodic");
        case BOUNDARY_REFLECTIVE:
            return PyUnicode_FromString("reflective");
        default:
            return PyUnicode_FromString("unknown");
    }
}

// 设置边界类型
static int Simulation_set_boundary(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete boundary attribute");
        return -1;
    }
    
    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "boundary must be a string");
        return -1;
    }
    
    const char *boundary_str = PyUnicode_AsUTF8(value);
    if (strcmp(boundary_str, "open") == 0) {
        self->sim->boundary = BOUNDARY_OPEN;
    } else if (strcmp(boundary_str, "periodic") == 0) {
        self->sim->boundary = BOUNDARY_PERIODIC;
    } else if (strcmp(boundary_str, "reflective") == 0) {
        self->sim->boundary = BOUNDARY_REFLECTIVE;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown boundary type");
        return -1;
    }
    
    return 0;
}

// 获取边界大小
static PyObject* Simulation_get_boundary_size(SimulationObject *self, void *closure) {
    return PyFloat_FromDouble(self->sim->boundary_size);
}

// 设置边界大小
static int Simulation_set_boundary_size(SimulationObject *self, PyObject *value, void *closure) {
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete boundary_size attribute");
        return -1;
    }
    
    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "boundary_size must be a float");
        return -1;
    }
    
    self->sim->boundary_size = PyFloat_AsDouble(value);
    return 0;
}

// 设置边界条件方法
static PyObject* Simulation_set_boundary_method(SimulationObject *self, PyObject *args, PyObject *kwargs) {
    char *boundary_str = "open";
    double size = 100.0;
    
    static char *kwlist[] = {"boundary", "size", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|sd", kwlist, &boundary_str, &size)) {
        return NULL;
    }
    
    BoundaryType boundary;
    if (strcmp(boundary_str, "open") == 0) {
        boundary = BOUNDARY_OPEN;
    } else if (strcmp(boundary_str, "periodic") == 0) {
        boundary = BOUNDARY_PERIODIC;
    } else if (strcmp(boundary_str, "reflective") == 0) {
        boundary = BOUNDARY_REFLECTIVE;
    } else {
        PyErr_SetString(PyExc_ValueError, "Unknown boundary type");
        return NULL;
    }
    
    sim_set_boundary(self->sim, boundary, size);
    Py_RETURN_NONE;
}

// 获取粒子列表
static PyObject* Simulation_get_particles(SimulationObject *self, void *closure) {
    PyObject *particle_list = PyList_New(self->sim->N);
    if (!particle_list) return NULL;
    
    for (int i = 0; i < self->sim->N; i++) {
        Particle *p = &self->sim->particles[i];
        
        // 创建粒子字典
        PyObject *particle_dict = PyDict_New();
        if (!particle_dict) {
            Py_DECREF(particle_list);
            return NULL;
        }
        
        PyDict_SetItemString(particle_dict, "m", PyFloat_FromDouble(p->m));
        PyDict_SetItemString(particle_dict, "x", PyFloat_FromDouble(p->x));
        PyDict_SetItemString(particle_dict, "y", PyFloat_FromDouble(p->y));
        PyDict_SetItemString(particle_dict, "z", PyFloat_FromDouble(p->z));
        PyDict_SetItemString(particle_dict, "vx", PyFloat_FromDouble(p->vx));
        PyDict_SetItemString(particle_dict, "vy", PyFloat_FromDouble(p->vy));
        PyDict_SetItemString(particle_dict, "vz", PyFloat_FromDouble(p->vz));
        PyDict_SetItemString(particle_dict, "r", PyFloat_FromDouble(p->r));
        PyDict_SetItemString(particle_dict, "id", PyLong_FromLong(p->id));
        
        PyList_SetItem(particle_list, i, particle_dict);
    }
    
    return particle_list;
}

// 移除粒子
static PyObject* Simulation_remove_particle(SimulationObject *self, PyObject *args) {
    int index;
    if (!PyArg_ParseTuple(args, "i", &index)) {
        return NULL;
    }
    
    if (index < 0 || index >= self->sim->N) {
        PyErr_SetString(PyExc_IndexError, "Particle index out of range");
        return NULL;
    }
    
    sim_remove_particle(self->sim, index);
    Py_RETURN_NONE;
}

// 获取单个粒子
static PyObject* Simulation_get_particle(SimulationObject *self, PyObject *args) {
    int index;
    if (!PyArg_ParseTuple(args, "i", &index)) {
        return NULL;
    }
    
    if (index < 0 || index >= self->sim->N) {
        PyErr_SetString(PyExc_IndexError, "Particle index out of range");
        return NULL;
    }
    
    Particle *p = &self->sim->particles[index];
    
    // 创建粒子字典
    PyObject *particle_dict = PyDict_New();
    if (!particle_dict) return NULL;
    
    PyDict_SetItemString(particle_dict, "m", PyFloat_FromDouble(p->m));
    PyDict_SetItemString(particle_dict, "x", PyFloat_FromDouble(p->x));
    PyDict_SetItemString(particle_dict, "y", PyFloat_FromDouble(p->y));
    PyDict_SetItemString(particle_dict, "z", PyFloat_FromDouble(p->z));
    PyDict_SetItemString(particle_dict, "vx", PyFloat_FromDouble(p->vx));
    PyDict_SetItemString(particle_dict, "vy", PyFloat_FromDouble(p->vy));
    PyDict_SetItemString(particle_dict, "vz", PyFloat_FromDouble(p->vz));
    PyDict_SetItemString(particle_dict, "r", PyFloat_FromDouble(p->r));
    PyDict_SetItemString(particle_dict, "id", PyLong_FromLong(p->id));
    
    return particle_dict;
}

// 计算动量
static PyObject* Simulation_calculate_momentum(SimulationObject *self, PyObject *args) {
    double px = 0.0, py = 0.0, pz = 0.0;
    
    for (int i = 0; i < self->sim->N; i++) {
        Particle *p = &self->sim->particles[i];
        px += p->m * p->vx;
        py += p->m * p->vy;
        pz += p->m * p->vz;
    }
    
    // 创建numpy数组
    npy_intp dims[1] = {3};
    PyObject *momentum_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!momentum_array) return NULL;
    
    double *data = (double*)PyArray_DATA((PyArrayObject*)momentum_array);
    data[0] = px;
    data[1] = py;
    data[2] = pz;
    
    return momentum_array;
}

// 计算能量
static PyObject* Simulation_calculate_energy(SimulationObject *self, PyObject *args) {
    double energy = calculate_energy(self->sim);
    return PyFloat_FromDouble(energy);
}

// 设置精确时间控制
static PyObject* Simulation_set_exact_finish_time(SimulationObject *self, PyObject *args) {
    int exact;
    if (!PyArg_ParseTuple(args, "i", &exact)) {
        return NULL;
    }
    
    sim_set_exact_finish_time(self->sim, exact);
    Py_RETURN_NONE;
}

// 方法定义
static PyMethodDef Simulation_methods[] = {
    {"add", (PyCFunction)Simulation_add_particle, METH_VARARGS | METH_KEYWORDS,
     "Add a particle to the simulation"},
    {"remove", (PyCFunction)Simulation_remove_particle, METH_VARARGS,
     "Remove a particle from the simulation"},
    {"get_particle", (PyCFunction)Simulation_get_particle, METH_VARARGS,
     "Get a single particle by index"},
    {"step", (PyCFunction)Simulation_step, METH_NOARGS,
     "Perform one integration step"},
    {"integrate", (PyCFunction)Simulation_integrate, METH_VARARGS,
     "Integrate to specified time"},
    {"calculate_energy", (PyCFunction)Simulation_calculate_energy, METH_NOARGS,
     "Calculate total energy of the system"},
    {"calculate_momentum", (PyCFunction)Simulation_calculate_momentum, METH_NOARGS,
     "Calculate total momentum of the system"},
    {"set_boundary", (PyCFunction)Simulation_set_boundary_method, METH_VARARGS | METH_KEYWORDS,
     "Set boundary conditions"},
     {"set_exact_finish_time", (PyCFunction)Simulation_set_exact_finish_time, METH_VARARGS,
     "Set exact finish time control"},
    {NULL}  /* Sentinel */
};

// 属性定义
static PyGetSetDef Simulation_getsetters[] = {
    {"N", (getter)Simulation_get_N, NULL, "Number of particles", NULL},
    {"t", (getter)Simulation_get_t, (setter)Simulation_set_t, "Current time", NULL},
    {"dt", (getter)Simulation_get_dt, (setter)Simulation_set_dt, "Time step", NULL},
    {"G", (getter)Simulation_get_G, (setter)Simulation_set_G, "Gravitational constant", NULL},
    {"softening", (getter)Simulation_get_softening, (setter)Simulation_set_softening, "Softening parameter", NULL},
    {"theta", (getter)Simulation_get_theta, (setter)Simulation_set_theta, "Barnes-Hut opening angle parameter", NULL},
    {"integrator", (getter)Simulation_get_integrator, (setter)Simulation_set_integrator, "Integrator type", NULL},
    {"collision", (getter)Simulation_get_collision, (setter)Simulation_set_collision, "Collision type", NULL},
    {"collision_detection", (getter)Simulation_get_collision_detection, (setter)Simulation_set_collision_detection, "Collision detection method", NULL},
    {"gravity_method", (getter)Simulation_get_gravity_method, (setter)Simulation_set_gravity_method, "Gravity calculation method", NULL},
    {"boundary", (getter)Simulation_get_boundary, (setter)Simulation_set_boundary, "Boundary type", NULL},
    {"boundary_size", (getter)Simulation_get_boundary_size, (setter)Simulation_set_boundary_size, "Boundary size", NULL},
    {"particles", (getter)Simulation_get_particles, NULL, "List of particles", NULL},
    {NULL}  /* Sentinel */
};

// 类型定义
static PyTypeObject SimulationType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "simple_rebound.core.Simulation",
    .tp_doc = "Simulation object",
    .tp_basicsize = sizeof(SimulationObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = Simulation_new,
    .tp_dealloc = (destructor)Simulation_dealloc,
    .tp_methods = Simulation_methods,
    .tp_getset = Simulation_getsetters,
};

// 模块方法
static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

// 模块定义
static struct PyModuleDef coremodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "core",
    .m_doc = "Simple Rebound core module",
    .m_size = -1,
    .m_methods = module_methods,
};

// 模块初始化
PyMODINIT_FUNC PyInit_core(void) {
    PyObject *m;
    
    if (PyType_Ready(&SimulationType) < 0)
        return NULL;
    
    m = PyModule_Create(&coremodule);
    if (m == NULL)
        return NULL;
    
    // 导入numpy
    import_array();
    
    Py_INCREF(&SimulationType);
    if (PyModule_AddObject(m, "Simulation", (PyObject *)&SimulationType) < 0) {
        Py_DECREF(&SimulationType);
        Py_DECREF(m);
        return NULL;
    }
    
    return m;
}
