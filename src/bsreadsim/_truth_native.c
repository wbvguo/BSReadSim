#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <math.h>
#include <string.h>


typedef struct {
    char *data;
    Py_ssize_t length;
    Py_ssize_t capacity;
} JsonBuffer;


typedef int (*JsonItemAppender)(JsonBuffer *, PyObject *);


#define APPEND_LITERAL(buffer, literal) \
    json_buffer_append((buffer), (literal), (Py_ssize_t)(sizeof(literal) - 1))


static void
json_buffer_free(JsonBuffer *buffer)
{
    PyMem_Free(buffer->data);
    buffer->data = NULL;
    buffer->length = 0;
    buffer->capacity = 0;
}


static int
json_buffer_reserve(JsonBuffer *buffer, Py_ssize_t additional)
{
    Py_ssize_t required;
    Py_ssize_t capacity;
    char *updated;

    if (additional < 0 || buffer->length > PY_SSIZE_T_MAX - additional) {
        PyErr_NoMemory();
        return -1;
    }
    required = buffer->length + additional;
    if (required <= buffer->capacity) {
        return 0;
    }

    capacity = buffer->capacity > 0 ? buffer->capacity : 4096;
    while (capacity < required) {
        if (capacity > PY_SSIZE_T_MAX / 2) {
            capacity = required;
            break;
        }
        capacity *= 2;
    }
    updated = PyMem_Realloc(buffer->data, (size_t)capacity);
    if (updated == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    buffer->data = updated;
    buffer->capacity = capacity;
    return 0;
}


static int
json_buffer_append(JsonBuffer *buffer, const char *value, Py_ssize_t length)
{
    if (length == 0) {
        return 0;
    }
    if (json_buffer_reserve(buffer, length) < 0) {
        return -1;
    }
    memcpy(buffer->data + buffer->length, value, (size_t)length);
    buffer->length += length;
    return 0;
}


static int
json_buffer_append_char(JsonBuffer *buffer, char value)
{
    return json_buffer_append(buffer, &value, 1);
}


static int
append_json_unicode(JsonBuffer *buffer, PyObject *value)
{
    const char *text;
    Py_ssize_t length;
    Py_ssize_t run_start;
    Py_ssize_t offset;

    if (!PyUnicode_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "truth JSON string field must be str");
        return -1;
    }
    text = PyUnicode_AsUTF8AndSize(value, &length);
    if (text == NULL) {
        return -1;
    }
    if (json_buffer_append_char(buffer, '"') < 0) {
        return -1;
    }

    run_start = 0;
    for (offset = 0; offset < length; ++offset) {
        const unsigned char byte = (unsigned char)text[offset];
        const char *escape = NULL;
        Py_ssize_t escape_length = 0;
        char unicode_escape[6];

        switch (byte) {
            case '"': escape = "\\\""; escape_length = 2; break;
            case '\\': escape = "\\\\"; escape_length = 2; break;
            case '\b': escape = "\\b"; escape_length = 2; break;
            case '\f': escape = "\\f"; escape_length = 2; break;
            case '\n': escape = "\\n"; escape_length = 2; break;
            case '\r': escape = "\\r"; escape_length = 2; break;
            case '\t': escape = "\\t"; escape_length = 2; break;
            default:
                if (byte < 0x20) {
                    static const char hex[] = "0123456789abcdef";
                    unicode_escape[0] = '\\';
                    unicode_escape[1] = 'u';
                    unicode_escape[2] = '0';
                    unicode_escape[3] = '0';
                    unicode_escape[4] = hex[(byte >> 4) & 0x0f];
                    unicode_escape[5] = hex[byte & 0x0f];
                    escape = unicode_escape;
                    escape_length = 6;
                }
                break;
        }
        if (escape == NULL) {
            continue;
        }
        if (json_buffer_append(buffer, text + run_start, offset - run_start) < 0
                || json_buffer_append(buffer, escape, escape_length) < 0) {
            return -1;
        }
        run_start = offset + 1;
    }
    if (json_buffer_append(buffer, text + run_start, length - run_start) < 0
            || json_buffer_append_char(buffer, '"') < 0) {
        return -1;
    }
    return 0;
}


static int
append_json_integer(JsonBuffer *buffer, PyObject *value)
{
    int overflow = 0;
    long long signed_value;
    char text[32];
    int length;

    if (!PyLong_Check(value) || value == Py_True || value == Py_False) {
        PyErr_SetString(PyExc_TypeError, "truth JSON integer field must be int");
        return -1;
    }
    signed_value = PyLong_AsLongLongAndOverflow(value, &overflow);
    if (signed_value == -1 && PyErr_Occurred()) {
        return -1;
    }
    if (overflow == 0) {
        length = PyOS_snprintf(text, sizeof(text), "%lld", signed_value);
    } else if (overflow > 0) {
        const unsigned long long unsigned_value =
            PyLong_AsUnsignedLongLong(value);
        if (unsigned_value == (unsigned long long)-1 && PyErr_Occurred()) {
            return -1;
        }
        length = PyOS_snprintf(text, sizeof(text), "%llu", unsigned_value);
    } else {
        PyErr_SetString(PyExc_OverflowError, "truth JSON integer is below int64");
        return -1;
    }
    if (length < 0 || length >= (int)sizeof(text)) {
        PyErr_SetString(PyExc_SystemError, "truth JSON integer formatting failed");
        return -1;
    }
    return json_buffer_append(buffer, text, (Py_ssize_t)length);
}


static int
append_json_unsigned_integer(JsonBuffer *buffer, unsigned long long value)
{
    char text[32];
    const int length = PyOS_snprintf(text, sizeof(text), "%llu", value);
    if (length < 0 || length >= (int)sizeof(text)) {
        PyErr_SetString(PyExc_SystemError, "truth JSON integer formatting failed");
        return -1;
    }
    return json_buffer_append(buffer, text, (Py_ssize_t)length);
}


static int
append_json_float(JsonBuffer *buffer, PyObject *value)
{
    double converted;
    char *text;
    int result;

    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "truth JSON probability field must be float");
        return -1;
    }
    converted = PyFloat_AsDouble(value);
    if (converted == -1.0 && PyErr_Occurred()) {
        return -1;
    }
    if (!isfinite(converted)) {
        PyErr_SetString(
            PyExc_ValueError,
            "Out of range float values are not JSON compliant"
        );
        return -1;
    }
    text = PyOS_double_to_string(
        converted,
        'r',
        0,
        Py_DTSF_ADD_DOT_0,
        NULL
    );
    if (text == NULL) {
        return -1;
    }
    result = json_buffer_append(buffer, text, (Py_ssize_t)strlen(text));
    PyMem_Free(text);
    return result;
}


static int
append_json_bool(JsonBuffer *buffer, PyObject *value)
{
    if (value == Py_True) {
        return APPEND_LITERAL(buffer, "true");
    }
    if (value == Py_False) {
        return APPEND_LITERAL(buffer, "false");
    }
    PyErr_SetString(PyExc_TypeError, "truth JSON boolean field must be bool");
    return -1;
}


static int
append_json_optional_bool(JsonBuffer *buffer, PyObject *value)
{
    if (value == Py_None) {
        return APPEND_LITERAL(buffer, "null");
    }
    return append_json_bool(buffer, value);
}


static int
append_json_optional_integer(JsonBuffer *buffer, PyObject *value)
{
    if (value == Py_None) {
        return APPEND_LITERAL(buffer, "null");
    }
    return append_json_integer(buffer, value);
}


static int
append_string_attribute(JsonBuffer *buffer, PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_unicode(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_integer_attribute(JsonBuffer *buffer, PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_integer(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_optional_integer_attribute(
    JsonBuffer *buffer,
    PyObject *object,
    const char *name
)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_optional_integer(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_float_attribute(JsonBuffer *buffer, PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_float(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_bool_attribute(JsonBuffer *buffer, PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_bool(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_optional_bool_attribute(
    JsonBuffer *buffer,
    PyObject *object,
    const char *name
)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_json_optional_bool(buffer, value);
    Py_DECREF(value);
    return result;
}


static int
append_enum_name_attribute(
    JsonBuffer *buffer,
    PyObject *object,
    const char *name
)
{
    PyObject *enumeration = PyObject_GetAttrString(object, name);
    PyObject *enum_name;
    int result;
    if (enumeration == NULL) {
        return -1;
    }
    enum_name = PyObject_GetAttrString(enumeration, "name");
    Py_DECREF(enumeration);
    if (enum_name == NULL) {
        return -1;
    }
    result = append_json_unicode(buffer, enum_name);
    Py_DECREF(enum_name);
    return result;
}


static int
append_base_string_attribute(
    JsonBuffer *buffer,
    PyObject *object,
    const char *name
)
{
    static const char bases[] = "ACGTN";
    PyObject *value = PyObject_GetAttrString(object, name);
    const unsigned char *data;
    Py_ssize_t length;
    Py_ssize_t offset;
    int result = -1;

    if (value == NULL) {
        return -1;
    }
    if (!PyBytes_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "variant bases must be bytes");
        goto done;
    }
    data = (const unsigned char *)PyBytes_AS_STRING(value);
    length = PyBytes_GET_SIZE(value);
    if (json_buffer_append_char(buffer, '"') < 0) {
        goto done;
    }
    for (offset = 0; offset < length; ++offset) {
        if (data[offset] > 4) {
            PyErr_SetString(
                PyExc_ValueError,
                "variant event bases are outside protocol encoding"
            );
            goto done;
        }
        if (json_buffer_append_char(buffer, bases[data[offset]]) < 0) {
            goto done;
        }
    }
    if (json_buffer_append_char(buffer, '"') < 0) {
        goto done;
    }
    result = 0;

done:
    Py_DECREF(value);
    return result;
}


static int
append_sequence_value(
    JsonBuffer *buffer,
    PyObject *value,
    JsonItemAppender append_item
)
{
    PyObject *sequence;
    Py_ssize_t index;
    Py_ssize_t length;
    int result = -1;

    sequence = PySequence_Fast(value, "truth JSON collection must be a sequence");
    if (sequence == NULL) {
        return -1;
    }
    if (json_buffer_append_char(buffer, '[') < 0) {
        goto done;
    }
    length = PySequence_Fast_GET_SIZE(sequence);
    for (index = 0; index < length; ++index) {
        if (index > 0 && json_buffer_append_char(buffer, ',') < 0) {
            goto done;
        }
        if (append_item(buffer, PySequence_Fast_GET_ITEM(sequence, index)) < 0) {
            goto done;
        }
    }
    if (json_buffer_append_char(buffer, ']') < 0) {
        goto done;
    }
    result = 0;

done:
    Py_DECREF(sequence);
    return result;
}


static int
append_sequence_attribute(
    JsonBuffer *buffer,
    PyObject *object,
    const char *name,
    JsonItemAppender append_item
)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    int result;
    if (value == NULL) {
        return -1;
    }
    result = append_sequence_value(buffer, value, append_item);
    Py_DECREF(value);
    return result;
}


typedef struct {
    PyObject *reference_positions;
    PyObject *variant_event_ids;
    PyObject *site_indices;
    PyObject *methylated;
    PyObject *oriented_bases;
    PyObject *post_conversion_bases;
    PyObject *final_bases;
    PyObject *attempted;
    PyObject *succeeded;
    PyObject *error_flags;
    PyObject *quality_phreds;
    Py_ssize_t length;
} CompactAnnotationColumns;


static void
compact_annotation_columns_clear(CompactAnnotationColumns *columns)
{
    Py_XDECREF(columns->reference_positions);
    Py_XDECREF(columns->variant_event_ids);
    Py_XDECREF(columns->site_indices);
    Py_XDECREF(columns->methylated);
    Py_XDECREF(columns->oriented_bases);
    Py_XDECREF(columns->post_conversion_bases);
    Py_XDECREF(columns->final_bases);
    Py_XDECREF(columns->attempted);
    Py_XDECREF(columns->succeeded);
    Py_XDECREF(columns->error_flags);
    Py_XDECREF(columns->quality_phreds);
    memset(columns, 0, sizeof(*columns));
}


static PyObject *
get_fast_sequence_attribute(PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    PyObject *sequence;
    if (value == NULL) {
        return NULL;
    }
    sequence = PySequence_Fast(value, "compact annotation column must be a sequence");
    Py_DECREF(value);
    return sequence;
}


static PyObject *
get_bytes_attribute(PyObject *object, const char *name)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    if (value == NULL) {
        return NULL;
    }
    if (!PyBytes_Check(value)) {
        Py_DECREF(value);
        PyErr_SetString(PyExc_TypeError, "compact annotation base column must be bytes");
        return NULL;
    }
    return value;
}


static int
load_compact_annotation_columns(
    PyObject *object,
    CompactAnnotationColumns *columns
)
{
    PyObject *sequence_columns[8];
    PyObject *byte_columns[3];
    Py_ssize_t index;

    memset(columns, 0, sizeof(*columns));
    columns->reference_positions = get_fast_sequence_attribute(
        object, "reference_positions"
    );
    columns->variant_event_ids = get_fast_sequence_attribute(
        object, "variant_event_ids"
    );
    columns->site_indices = get_fast_sequence_attribute(object, "site_indices");
    columns->methylated = get_fast_sequence_attribute(object, "methylated");
    columns->oriented_bases = get_bytes_attribute(object, "oriented_bases");
    columns->post_conversion_bases = get_bytes_attribute(
        object, "post_conversion_bases"
    );
    columns->final_bases = get_bytes_attribute(object, "final_bases");
    columns->attempted = get_fast_sequence_attribute(object, "attempted");
    columns->succeeded = get_fast_sequence_attribute(object, "succeeded");
    columns->error_flags = get_fast_sequence_attribute(object, "error_flags");
    columns->quality_phreds = get_fast_sequence_attribute(
        object, "quality_phreds"
    );
    if (columns->reference_positions == NULL
            || columns->variant_event_ids == NULL
            || columns->site_indices == NULL
            || columns->methylated == NULL
            || columns->oriented_bases == NULL
            || columns->post_conversion_bases == NULL
            || columns->final_bases == NULL
            || columns->attempted == NULL
            || columns->succeeded == NULL
            || columns->error_flags == NULL
            || columns->quality_phreds == NULL) {
        compact_annotation_columns_clear(columns);
        return -1;
    }

    columns->length = PyBytes_GET_SIZE(columns->final_bases);
    sequence_columns[0] = columns->reference_positions;
    sequence_columns[1] = columns->variant_event_ids;
    sequence_columns[2] = columns->site_indices;
    sequence_columns[3] = columns->methylated;
    sequence_columns[4] = columns->attempted;
    sequence_columns[5] = columns->succeeded;
    sequence_columns[6] = columns->error_flags;
    sequence_columns[7] = columns->quality_phreds;
    for (index = 0; index < 8; ++index) {
        if (PySequence_Fast_GET_SIZE(sequence_columns[index]) != columns->length) {
            compact_annotation_columns_clear(columns);
            PyErr_SetString(
                PyExc_ValueError,
                "compact annotation columns have different lengths"
            );
            return -1;
        }
    }
    byte_columns[0] = columns->oriented_bases;
    byte_columns[1] = columns->post_conversion_bases;
    byte_columns[2] = columns->final_bases;
    for (index = 0; index < 3; ++index) {
        if (PyBytes_GET_SIZE(byte_columns[index]) != columns->length) {
            compact_annotation_columns_clear(columns);
            PyErr_SetString(
                PyExc_ValueError,
                "compact annotation columns have different lengths"
            );
            return -1;
        }
    }
    return 0;
}


static int append_annotation(JsonBuffer *buffer, PyObject *annotation);
static int append_mate(JsonBuffer *buffer, PyObject *mate);
static int append_site_state(JsonBuffer *buffer, PyObject *site);
static int append_variant_event(JsonBuffer *buffer, PyObject *event);


static int
append_compact_annotation(
    JsonBuffer *buffer,
    const CompactAnnotationColumns *columns,
    Py_ssize_t offset
)
{
    const unsigned char oriented_base = (unsigned char)
        PyBytes_AS_STRING(columns->oriented_bases)[offset];
    const unsigned char post_conversion_base = (unsigned char)
        PyBytes_AS_STRING(columns->post_conversion_bases)[offset];
    const unsigned char final_base = (unsigned char)
        PyBytes_AS_STRING(columns->final_bases)[offset];

    if (APPEND_LITERAL(buffer, "{\"conversion_attempted\":") < 0
            || append_json_bool(
                buffer,
                PySequence_Fast_GET_ITEM(columns->attempted, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"conversion_succeeded\":") < 0
            || append_json_bool(
                buffer,
                PySequence_Fast_GET_ITEM(columns->succeeded, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"final_base\":") < 0
            || append_json_unsigned_integer(buffer, final_base) < 0
            || APPEND_LITERAL(buffer, ",\"methylated\":") < 0
            || append_json_optional_bool(
                buffer,
                PySequence_Fast_GET_ITEM(columns->methylated, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"oriented_base\":") < 0
            || append_json_unsigned_integer(buffer, oriented_base) < 0
            || APPEND_LITERAL(buffer, ",\"post_conversion_base\":") < 0
            || append_json_unsigned_integer(buffer, post_conversion_base) < 0
            || APPEND_LITERAL(buffer, ",\"quality_phred\":") < 0
            || append_json_integer(
                buffer,
                PySequence_Fast_GET_ITEM(columns->quality_phreds, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"read_offset\":") < 0
            || append_json_unsigned_integer(
                buffer, (unsigned long long)offset
            ) < 0
            || APPEND_LITERAL(buffer, ",\"reference_pos\":") < 0
            || append_json_integer(
                buffer,
                PySequence_Fast_GET_ITEM(columns->reference_positions, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"sequencing_error\":") < 0
            || append_json_bool(
                buffer,
                PySequence_Fast_GET_ITEM(columns->error_flags, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"site_index\":") < 0
            || append_json_optional_integer(
                buffer,
                PySequence_Fast_GET_ITEM(columns->site_indices, offset)
            ) < 0
            || APPEND_LITERAL(buffer, ",\"variant_event_id\":") < 0
            || append_json_integer(
                buffer,
                PySequence_Fast_GET_ITEM(columns->variant_event_ids, offset)
            ) < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


static int
append_annotations_attribute(JsonBuffer *buffer, PyObject *mate)
{
    PyObject *annotations = PyObject_GetAttrString(mate, "annotations");
    CompactAnnotationColumns columns;
    Py_ssize_t offset;
    int result = -1;

    if (annotations == NULL) {
        return -1;
    }
    if (PyTuple_Check(annotations) || PyList_Check(annotations)) {
        result = append_sequence_value(buffer, annotations, append_annotation);
        Py_DECREF(annotations);
        return result;
    }
    if (load_compact_annotation_columns(annotations, &columns) < 0) {
        Py_DECREF(annotations);
        return -1;
    }
    Py_DECREF(annotations);

    if (json_buffer_append_char(buffer, '[') < 0) {
        goto done;
    }
    for (offset = 0; offset < columns.length; ++offset) {
        if (offset > 0 && json_buffer_append_char(buffer, ',') < 0) {
            goto done;
        }
        if (append_compact_annotation(buffer, &columns, offset) < 0) {
            goto done;
        }
    }
    if (json_buffer_append_char(buffer, ']') < 0) {
        goto done;
    }
    result = 0;

done:
    compact_annotation_columns_clear(&columns);
    return result;
}


static int
append_annotation(JsonBuffer *buffer, PyObject *annotation)
{
    if (APPEND_LITERAL(buffer, "{\"conversion_attempted\":") < 0
            || append_bool_attribute(
                buffer, annotation, "conversion_attempted"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"conversion_succeeded\":") < 0
            || append_bool_attribute(
                buffer, annotation, "conversion_succeeded"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"final_base\":") < 0
            || append_integer_attribute(buffer, annotation, "final_base") < 0
            || APPEND_LITERAL(buffer, ",\"methylated\":") < 0
            || append_optional_bool_attribute(
                buffer, annotation, "methylated"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"oriented_base\":") < 0
            || append_integer_attribute(buffer, annotation, "oriented_base") < 0
            || APPEND_LITERAL(buffer, ",\"post_conversion_base\":") < 0
            || append_integer_attribute(
                buffer, annotation, "post_conversion_base"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"quality_phred\":") < 0
            || append_integer_attribute(buffer, annotation, "quality_phred") < 0
            || APPEND_LITERAL(buffer, ",\"read_offset\":") < 0
            || append_integer_attribute(buffer, annotation, "read_offset") < 0
            || APPEND_LITERAL(buffer, ",\"reference_pos\":") < 0
            || append_integer_attribute(buffer, annotation, "reference_pos") < 0
            || APPEND_LITERAL(buffer, ",\"sequencing_error\":") < 0
            || append_bool_attribute(buffer, annotation, "sequencing_error") < 0
            || APPEND_LITERAL(buffer, ",\"site_index\":") < 0
            || append_optional_integer_attribute(
                buffer, annotation, "site_index"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"variant_event_id\":") < 0
            || append_integer_attribute(
                buffer, annotation, "variant_event_id"
            ) < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


static int
append_mate(JsonBuffer *buffer, PyObject *mate)
{
    if (APPEND_LITERAL(buffer, "{\"annotations\":") < 0
            || append_annotations_attribute(buffer, mate) < 0
            || APPEND_LITERAL(buffer, ",\"conversion_mode\":") < 0
            || append_enum_name_attribute(buffer, mate, "conversion_mode") < 0
            || APPEND_LITERAL(buffer, ",\"mate_index\":") < 0
            || append_integer_attribute(buffer, mate, "mate_index") < 0
            || APPEND_LITERAL(buffer, ",\"reference_end\":") < 0
            || append_integer_attribute(buffer, mate, "reference_end") < 0
            || APPEND_LITERAL(buffer, ",\"reference_start\":") < 0
            || append_integer_attribute(buffer, mate, "reference_start") < 0
            || APPEND_LITERAL(buffer, ",\"reverse_complement\":") < 0
            || append_bool_attribute(buffer, mate, "reverse_complement") < 0
            || APPEND_LITERAL(buffer, ",\"sequence\":") < 0
            || append_string_attribute(buffer, mate, "sequence") < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


static int
append_site_state(JsonBuffer *buffer, PyObject *site)
{
    if (APPEND_LITERAL(buffer, "{\"allele\":") < 0
            || append_enum_name_attribute(buffer, site, "allele") < 0
            || APPEND_LITERAL(buffer, ",\"context\":") < 0
            || append_enum_name_attribute(buffer, site, "context") < 0
            || APPEND_LITERAL(buffer, ",\"methylated\":") < 0
            || append_bool_attribute(buffer, site, "methylated") < 0
            || APPEND_LITERAL(buffer, ",\"probability\":") < 0
            || append_float_attribute(buffer, site, "probability") < 0
            || APPEND_LITERAL(buffer, ",\"reference_pos\":") < 0
            || append_integer_attribute(buffer, site, "reference_pos") < 0
            || APPEND_LITERAL(buffer, ",\"site_index\":") < 0
            || append_integer_attribute(buffer, site, "site_index") < 0
            || APPEND_LITERAL(buffer, ",\"source\":") < 0
            || append_enum_name_attribute(buffer, site, "source") < 0
            || APPEND_LITERAL(buffer, ",\"template_offset\":") < 0
            || append_integer_attribute(buffer, site, "template_offset") < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


static int
append_variant_event(JsonBuffer *buffer, PyObject *event)
{
    if (APPEND_LITERAL(buffer, "{\"alt_bases\":") < 0
            || append_base_string_attribute(buffer, event, "alt_bases") < 0
            || APPEND_LITERAL(buffer, ",\"event_id\":") < 0
            || append_integer_attribute(buffer, event, "event_id") < 0
            || APPEND_LITERAL(buffer, ",\"kind\":") < 0
            || append_enum_name_attribute(buffer, event, "kind") < 0
            || APPEND_LITERAL(buffer, ",\"phased_haplotype\":") < 0
            || append_integer_attribute(buffer, event, "phased_haplotype") < 0
            || APPEND_LITERAL(buffer, ",\"ref_bases\":") < 0
            || append_base_string_attribute(buffer, event, "ref_bases") < 0
            || APPEND_LITERAL(buffer, ",\"reference_end\":") < 0
            || append_integer_attribute(buffer, event, "reference_end") < 0
            || APPEND_LITERAL(buffer, ",\"reference_start\":") < 0
            || append_integer_attribute(buffer, event, "reference_start") < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


static int
append_fragment(JsonBuffer *buffer, PyObject *fragment)
{
    if (APPEND_LITERAL(buffer, "{\"contig\":") < 0
            || append_string_attribute(buffer, fragment, "contig_name") < 0
            || APPEND_LITERAL(buffer, ",\"fragment_conversion_mode\":") < 0
            || append_enum_name_attribute(
                buffer, fragment, "fragment_conversion_mode"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"fragment_ordinal\":") < 0
            || append_integer_attribute(
                buffer, fragment, "fragment_ordinal"
            ) < 0
            || APPEND_LITERAL(buffer, ",\"haplotype\":") < 0
            || append_integer_attribute(buffer, fragment, "haplotype") < 0
            || APPEND_LITERAL(buffer, ",\"mates\":") < 0
            || append_sequence_attribute(buffer, fragment, "mates", append_mate) < 0
            || APPEND_LITERAL(buffer, ",\"site_states\":") < 0
            || append_sequence_attribute(
                buffer, fragment, "site_states", append_site_state
            ) < 0
            || APPEND_LITERAL(buffer, ",\"variant_events\":") < 0
            || append_sequence_attribute(
                buffer, fragment, "variant_events", append_variant_event
            ) < 0
            || json_buffer_append_char(buffer, '}') < 0) {
        return -1;
    }
    return 0;
}


PyObject *
bsreadsim_native_truth_json_bytes(PyObject *self, PyObject *args)
{
    PyObject *fragment;
    int newline = 0;
    JsonBuffer buffer = {NULL, 0, 0};
    PyObject *result;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "O|p:canonical_truth_json_bytes",
            &fragment,
            &newline
        )) {
        return NULL;
    }
    if (append_fragment(&buffer, fragment) < 0) {
        json_buffer_free(&buffer);
        return NULL;
    }
    if (newline && json_buffer_append_char(&buffer, '\n') < 0) {
        json_buffer_free(&buffer);
        return NULL;
    }
    result = PyBytes_FromStringAndSize(buffer.data, buffer.length);
    json_buffer_free(&buffer);
    return result;
}
