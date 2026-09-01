#ifndef BSREADSIM_CEXT_API_H
#define BSREADSIM_CEXT_API_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdint.h>


extern PyObject *fragment_type;
extern PyObject *variant_type;
extern PyObject *methylation_site_type;
extern PyObject *site_reference_type;
extern PyObject *mate_type;
extern PyObject *capture_strands[3];
extern PyObject *variant_sources[4];
extern PyObject *variant_kinds[4];
extern PyObject *methylation_contexts[16];
extern PyObject *methylation_sources[5];
extern PyObject *methylation_alleles[3];

int initialize_protocol_types(void);
int tuple_set_unsigned(PyObject *tuple, Py_ssize_t index, uint64_t value);
int tuple_set_signed(PyObject *tuple, Py_ssize_t index, int64_t value);
int tuple_set_borrowed(PyObject *tuple, Py_ssize_t index, PyObject *value);
PyObject *call_record(PyObject *record_type, PyObject *arguments);

PyObject *bsreadsim_cext_format_sam_batch(PyObject *self, PyObject *args);
PyObject *bsreadsim_cext_format_sam_columns(PyObject *self, PyObject *args);
PyObject *bsreadsim_cext_validate_protocol_batch_columns(
    PyObject *self,
    PyObject *args
);
PyObject *bsreadsim_cext_pack_protocol_common_columns(
    PyObject *self,
    PyObject *args
);
PyObject *bsreadsim_cext_decode_protocol_fragments(
    PyObject *self,
    PyObject *args
);

PyMODINIT_FUNC PyInit__cext(void);

#endif
