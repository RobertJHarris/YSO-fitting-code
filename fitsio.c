 /* test obtained from http://blog.numerix-dsp.com/2013/08/how-to-pass-c-array-to-python-solution.html */

#include <Python/Python.h>
#include <stdio.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include/numpy/arrayobject.h"

/* need to put all initialization / finalization stuff here because
 * it can only be called ONCE. otherwise it'll seg fault for mysterious 
 * reasons 
 */


void quit_python(){
  Py_Finalize();     /*for whatever reason, this causes a segfault....*/;
}

void init_numpy(){
  import_array(); /* import_array() is a macro that has a return statement in it. can't use it in main */
}

void initialize_python(){
  Py_Initialize ();  
  init_numpy();
}

void  *set_path(){
  /* set the system path */
  PyObject *sys_path; 
  PyObject *path; 
  sys_path = PySys_GetObject("path"); 
  if (sys_path == NULL) 
    return NULL; 
  path = PyString_FromString(".");
    if (path == NULL) 
      return NULL; 
  if (PyList_Append(sys_path, path) < 0) 
    return NULL; 
      Py_RETURN_NONE; 
}

void make_fits_model_im(char *filename, double *im, int *size, int dim, double pix_size)
{
  PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pValue;
  PyObject *py_filename;
  npy_intp dims[dim];
  PyObject *py_array;
  PyObject *py_pix;
  PyObject *path;
  PyObject *py_refpix;
  int i;

  for(i=0; i<dim;i++)
    dims[i] = size[i];
 
  if(set_path() == NULL){
    fprintf(stderr,"ERROR setting Python path\n");
    exit(-1);
  }

  pName = PyUnicode_FromString ("fitsio");
  // PyUnicode_FromString error checking here
  pModule = PyImport_Import (pName);                  // Load the module object
  if(pModule == NULL){
    fprintf(stderr,"ERROR importing module\n");
    exit(-1);
  } 

  // PyImport_Import error checking here
  Py_DECREF(pName);

  pFunc = PyObject_GetAttrString (pModule, "make_fits_model");        // pFunc is a new reference
  // PyObject_GetAttrString error checking here
  if(pModule == NULL){
    fprintf(stderr,"ERROR importing function\n");
    exit(-1);
  } 

  


  // Required for the C-API : http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
  //  fprintf(stderr,"[In fitsio.c] ]dim = %d \n",dim);
  //  fprintf(stderr,"[In fitsio.c] ]refpix = %d \n",size[0]/2+1);
  
  py_filename = PyUnicode_FromString (filename);
  py_array = PyArray_SimpleNewFromData (dim, dims, NPY_DOUBLE, im);
  py_pix   = PyFloat_FromDouble(pix_size);
  py_refpix   = PyFloat_FromDouble((size[1] % 2 == 0) ? size[1]/2+1 : (size[1]-1)/2+1);
  

  // PyArray_SimpleNewFromData error checking here

  //  return;
  pDict = PyModule_GetDict (pModule);                 // pDict is a borrowed reference
  // fprintf(stderr,"PyFloat %g\n",PyFloat_AS_DOUBLE(py_pix));
  pArgs = PyTuple_New (4);

  PyTuple_SetItem (pArgs, 0, py_filename);
  PyTuple_SetItem (pArgs, 1, py_array);
  PyTuple_SetItem (pArgs, 2, py_pix);
  PyTuple_SetItem (pArgs, 3, py_refpix);

  pFunc = PyDict_GetItemString (pDict, "make_fits_model"); // pFunc is also a borrowed reference

  if (PyCallable_Check (pFunc)) 
    {
      PyObject_CallObject (pFunc, pArgs);
    } else 
    {
      printf ("Function not callable !\n");
    }

  Py_DECREF (py_filename);
  Py_DECREF (py_array);                               // Clean up
  Py_DECREF (pModule);
  //  Py_DECREF (pDict);
  Py_DECREF (pFunc);

  return;
}
