#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdint.h>
#include <string.h>


#define MAX_CIGAR_OPERATION_LENGTH UINT64_C(268435455)


typedef struct {
    char *data;
    Py_ssize_t length;
    Py_ssize_t capacity;
} SamBuffer;


typedef struct {
    PyObject *sequence_object;
    PyObject *quality_object;
    const char *sequence;
    const char *quality;
    Py_ssize_t query_length;
    int quality_byte;
    long long mate_index;
    long long reference_start;
    long long reference_end;
    long long position;
    int reverse_complement;
    const char *base_state_codes;
    const unsigned char *read_summary;
    const unsigned char *fragment_summary;
    const char *fragment_realization;
    Py_ssize_t fragment_realization_length;
    SamBuffer cigar;
} MateAlignment;


static void
sam_buffer_clear(SamBuffer *buffer)
{
    PyMem_Free(buffer->data);
    buffer->data = NULL;
    buffer->length = 0;
    buffer->capacity = 0;
}


static int
sam_buffer_reserve(SamBuffer *buffer, Py_ssize_t additional)
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
sam_buffer_append(SamBuffer *buffer, const char *value, Py_ssize_t length)
{
    if (length == 0) {
        return 0;
    }
    if (sam_buffer_reserve(buffer, length) < 0) {
        return -1;
    }
    memcpy(buffer->data + buffer->length, value, (size_t)length);
    buffer->length += length;
    return 0;
}


static int
sam_buffer_append_char(SamBuffer *buffer, char value)
{
    return sam_buffer_append(buffer, &value, 1);
}


static int
sam_buffer_append_literal(SamBuffer *buffer, const char *value)
{
    return sam_buffer_append(buffer, value, (Py_ssize_t)strlen(value));
}


static int
sam_buffer_append_signed(SamBuffer *buffer, long long value)
{
    char text[32];
    const int length = PyOS_snprintf(text, sizeof(text), "%lld", value);
    if (length < 0 || length >= (int)sizeof(text)) {
        PyErr_SetString(PyExc_SystemError, "details SAM integer formatting failed");
        return -1;
    }
    return sam_buffer_append(buffer, text, (Py_ssize_t)length);
}


static int
sam_buffer_append_unsigned(SamBuffer *buffer, unsigned long long value)
{
    char text[32];
    const int length = PyOS_snprintf(text, sizeof(text), "%llu", value);
    if (length < 0 || length >= (int)sizeof(text)) {
        PyErr_SetString(PyExc_SystemError, "details SAM integer formatting failed");
        return -1;
    }
    return sam_buffer_append(buffer, text, (Py_ssize_t)length);
}


static int
sam_buffer_append_hex_unsigned(SamBuffer *buffer, unsigned long long value)
{
    char text[32];
    const int length = PyOS_snprintf(text, sizeof(text), "%llx", value);
    if (length < 0 || length >= (int)sizeof(text)) {
        PyErr_SetString(PyExc_SystemError, "details SAM hexadecimal formatting failed");
        return -1;
    }
    return sam_buffer_append(buffer, text, (Py_ssize_t)length);
}


static int
sam_buffer_append_unicode(SamBuffer *buffer, PyObject *value, const char *field)
{
    const char *text;
    Py_ssize_t length;
    if (!PyUnicode_Check(value)) {
        PyErr_Format(PyExc_TypeError, "details SAM %s must be str", field);
        return -1;
    }
    text = PyUnicode_AsUTF8AndSize(value, &length);
    if (text == NULL) {
        return -1;
    }
    return sam_buffer_append(buffer, text, length);
}


static int
get_integer_attribute(PyObject *object, const char *name, long long *result)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    long long converted;
    if (value == NULL) {
        return -1;
    }
    if (!PyLong_Check(value) || value == Py_True || value == Py_False) {
        PyErr_Format(PyExc_TypeError, "details SAM %s must be int", name);
        Py_DECREF(value);
        return -1;
    }
    converted = PyLong_AsLongLong(value);
    Py_DECREF(value);
    if (converted == -1 && PyErr_Occurred()) {
        return -1;
    }
    *result = converted;
    return 0;
}


static int
get_boolean_attribute(PyObject *object, const char *name, int *result)
{
    PyObject *value = PyObject_GetAttrString(object, name);
    if (value == NULL) {
        return -1;
    }
    if (value != Py_True && value != Py_False) {
        PyErr_Format(PyExc_TypeError, "details SAM %s must be bool", name);
        Py_DECREF(value);
        return -1;
    }
    *result = value == Py_True;
    Py_DECREF(value);
    return 0;
}


static int
read_reference_position(
    PyObject *positions,
    int direct,
    Py_ssize_t index,
    long long *result
)
{
    PyObject *item = PySequence_Fast_GET_ITEM(positions, index);
    PyObject *value = item;
    long long converted;
    if (!direct) {
        value = PyObject_GetAttrString(item, "reference_pos");
        if (value == NULL) {
            return -1;
        }
    }
    if (!PyLong_Check(value) || value == Py_True || value == Py_False) {
        PyErr_SetString(PyExc_TypeError, "details SAM reference positions must be integers");
        if (!direct) {
            Py_DECREF(value);
        }
        return -1;
    }
    converted = PyLong_AsLongLong(value);
    if (!direct) {
        Py_DECREF(value);
    }
    if (converted == -1 && PyErr_Occurred()) {
        return -1;
    }
    *result = converted;
    return 0;
}


static int
flush_cigar_operation(SamBuffer *buffer, char *operation, uint64_t *length)
{
    if (*length == 0) {
        return 0;
    }
    if (sam_buffer_append_unsigned(buffer, *length) < 0
            || sam_buffer_append_char(buffer, *operation) < 0) {
        return -1;
    }
    *operation = '\0';
    *length = 0;
    return 0;
}


static int
push_cigar_operation(
    SamBuffer *buffer,
    char *pending_operation,
    uint64_t *pending_length,
    char operation,
    uint64_t length
)
{
    while (length > 0) {
        uint64_t consumed;
        if (*pending_length > 0 && *pending_operation != operation) {
            if (flush_cigar_operation(
                    buffer, pending_operation, pending_length
                ) < 0) {
                return -1;
            }
        }
        if (*pending_length == 0) {
            *pending_operation = operation;
        }
        consumed = MAX_CIGAR_OPERATION_LENGTH - *pending_length;
        if (consumed > length) {
            consumed = length;
        }
        *pending_length += consumed;
        length -= consumed;
        if (*pending_length == MAX_CIGAR_OPERATION_LENGTH
                && flush_cigar_operation(
                    buffer, pending_operation, pending_length
                ) < 0) {
            return -1;
        }
    }
    return 0;
}


static void
mate_alignment_clear(MateAlignment *alignment)
{
    Py_XDECREF(alignment->sequence_object);
    Py_XDECREF(alignment->quality_object);
    alignment->sequence_object = NULL;
    alignment->quality_object = NULL;
    sam_buffer_clear(&alignment->cigar);
}


static int
build_mate_alignment(
    PyObject *mate,
    long long contig_length,
    MateAlignment *alignment
)
{
    PyObject *base_states = NULL;
    PyObject *position_values = NULL;
    PyObject *positions = NULL;
    Py_ssize_t sequence_length;
    Py_ssize_t quality_length;
    Py_ssize_t position_count;
    Py_ssize_t offset;
    int direct_positions = 0;
    int have_mapped = 0;
    long long previous_mapped = 0;
    long long mapped_start = 0;
    long long mapped_end = 0;
    long long declared_start;
    long long declared_end;
    uint64_t query_consumed = 0;
    char pending_operation = '\0';
    uint64_t pending_length = 0;
    int result = -1;

    memset(alignment, 0, sizeof(*alignment));
    if (get_integer_attribute(mate, "mate_index", &alignment->mate_index) < 0
            || get_boolean_attribute(
                mate, "reverse_complement", &alignment->reverse_complement
            ) < 0
            || get_integer_attribute(mate, "reference_start", &declared_start) < 0
            || get_integer_attribute(mate, "reference_end", &declared_end) < 0) {
        goto done;
    }
    alignment->sequence_object = PyObject_GetAttrString(mate, "sequence");
    alignment->quality_object = PyObject_GetAttrString(mate, "quality");
    if (alignment->sequence_object == NULL || alignment->quality_object == NULL) {
        goto done;
    }
    if (!PyUnicode_Check(alignment->sequence_object)
            || !PyUnicode_Check(alignment->quality_object)) {
        PyErr_SetString(PyExc_TypeError, "details SAM sequence and quality must be str");
        goto done;
    }
    alignment->sequence = PyUnicode_AsUTF8AndSize(
        alignment->sequence_object, &sequence_length
    );
    alignment->quality = PyUnicode_AsUTF8AndSize(
        alignment->quality_object, &quality_length
    );
    if (alignment->sequence == NULL || alignment->quality == NULL) {
        goto done;
    }
    if (sequence_length <= 0 || sequence_length != quality_length) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM sequence and quality lengths disagree"
        );
        goto done;
    }
    for (offset = 0; offset < sequence_length; ++offset) {
        const char base = alignment->sequence[offset];
        const unsigned char quality = (unsigned char)alignment->quality[offset];
        if (strchr("ACGTN", base) == NULL) {
            PyErr_SetString(PyExc_ValueError, "details SAM sequence is not A/C/G/T/N");
            goto done;
        }
        if (quality < 33 || quality > 126) {
            PyErr_SetString(PyExc_ValueError, "details SAM quality is not printable");
            goto done;
        }
    }
    alignment->query_length = sequence_length;

    base_states = PyObject_GetAttrString(mate, "base_states");
    if (base_states == NULL) {
        goto done;
    }
    position_values = PyObject_GetAttrString(base_states, "reference_positions");
    if (position_values != NULL) {
        direct_positions = 1;
        positions = PySequence_Fast(
            position_values,
            "details SAM compact reference positions must be a sequence"
        );
    } else {
        if (!PyErr_ExceptionMatches(PyExc_AttributeError)) {
            goto done;
        }
        PyErr_Clear();
        positions = PySequence_Fast(
            base_states,
            "details SAM base_states must be a sequence"
        );
    }
    if (positions == NULL) {
        goto done;
    }
    position_count = PySequence_Fast_GET_SIZE(positions);
    if (position_count != sequence_length) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM base_states must cover every read base"
        );
        goto done;
    }

    for (offset = 0; offset < position_count; ++offset) {
        const Py_ssize_t index = alignment->reverse_complement
            ? position_count - offset - 1
            : offset;
        long long position;
        if (read_reference_position(
                positions, direct_positions, index, &position
            ) < 0) {
            goto done;
        }
        if (position == -1) {
            if (push_cigar_operation(
                    &alignment->cigar,
                    &pending_operation,
                    &pending_length,
                    'I',
                    1
                ) < 0) {
                goto done;
            }
            query_consumed += 1;
            continue;
        }
        if (position < 0 || position >= contig_length) {
            PyErr_SetString(
                PyExc_ValueError,
                "details SAM reference position is outside its contig"
            );
            goto done;
        }
        if (have_mapped) {
            const long long gap = position - previous_mapped - 1;
            if (gap < 0) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "details SAM reference projection is not strictly increasing"
                );
                goto done;
            }
            if (gap > 0 && push_cigar_operation(
                    &alignment->cigar,
                    &pending_operation,
                    &pending_length,
                    'D',
                    (uint64_t)gap
                ) < 0) {
                goto done;
            }
        } else {
            mapped_start = position;
            have_mapped = 1;
        }
        if (push_cigar_operation(
                &alignment->cigar,
                &pending_operation,
                &pending_length,
                'M',
                1
            ) < 0) {
            goto done;
        }
        query_consumed += 1;
        previous_mapped = position;
        mapped_end = position + 1;
    }
    if (flush_cigar_operation(
            &alignment->cigar, &pending_operation, &pending_length
        ) < 0) {
        goto done;
    }
    if (query_consumed != (uint64_t)sequence_length) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM CIGAR does not consume the complete read"
        );
        goto done;
    }
    if (have_mapped) {
        if (declared_start != mapped_start || declared_end != mapped_end) {
            PyErr_SetString(
                PyExc_ValueError,
                "details SAM mate bounds disagree with its reference projection"
            );
            goto done;
        }
        alignment->reference_start = mapped_start;
        alignment->reference_end = mapped_end;
    } else {
        if (declared_start != declared_end
                || declared_start < 0
                || declared_start > contig_length) {
            PyErr_SetString(PyExc_ValueError, "details SAM insertion anchor is invalid");
            goto done;
        }
        alignment->reference_start = declared_start < contig_length
            ? declared_start
            : contig_length - 1;
        alignment->reference_end = alignment->reference_start + 1;
    }
    alignment->position = alignment->reference_start + 1;
    result = 0;

done:
    Py_XDECREF(positions);
    Py_XDECREF(position_values);
    Py_XDECREF(base_states);
    if (result < 0) {
        mate_alignment_clear(alignment);
    }
    return result;
}


static int
append_reference_forward_sequence(SamBuffer *buffer, const MateAlignment *alignment)
{
    Py_ssize_t offset;
    if (!alignment->reverse_complement) {
        return sam_buffer_append(
            buffer, alignment->sequence, alignment->query_length
        );
    }
    for (offset = alignment->query_length; offset > 0; --offset) {
        char base;
        switch (alignment->sequence[offset - 1]) {
            case 'A': base = 'T'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'T': base = 'A'; break;
            default: base = 'N'; break;
        }
        if (sam_buffer_append_char(buffer, base) < 0) {
            return -1;
        }
    }
    return 0;
}


static int
append_reference_forward_quality(SamBuffer *buffer, const MateAlignment *alignment)
{
    Py_ssize_t offset;
    if (alignment->quality == NULL) {
        if (sam_buffer_reserve(buffer, alignment->query_length) < 0) {
            return -1;
        }
        memset(
            buffer->data + buffer->length,
            alignment->quality_byte,
            (size_t)alignment->query_length
        );
        buffer->length += alignment->query_length;
        return 0;
    }
    if (!alignment->reverse_complement) {
        return sam_buffer_append(
            buffer, alignment->quality, alignment->query_length
        );
    }
    for (offset = alignment->query_length; offset > 0; --offset) {
        if (sam_buffer_append_char(buffer, alignment->quality[offset - 1]) < 0) {
            return -1;
        }
    }
    return 0;
}


static int
append_query_name(
    SamBuffer *buffer,
    PyObject *contig,
    long long reference_start,
    long long reference_end,
    unsigned long long fragment_ordinal
)
{
    const long long right = reference_end > reference_start
        ? reference_end
        : reference_start + 1;
    return sam_buffer_append_unicode(buffer, contig, "contig name") < 0
        || sam_buffer_append_char(buffer, ':') < 0
        || sam_buffer_append_signed(buffer, reference_start + 1) < 0
        || sam_buffer_append_char(buffer, '-') < 0
        || sam_buffer_append_signed(buffer, right) < 0
        || sam_buffer_append_char(buffer, ':') < 0
        || sam_buffer_append_hex_unsigned(buffer, fragment_ordinal) < 0
        ? -1
        : 0;
}


static uint16_t
load_u16_le_bytes(const unsigned char *value)
{
    return (uint16_t)value[0] | ((uint16_t)value[1] << 8);
}


static int
append_u16_summary_values(SamBuffer *buffer, const unsigned char *summary)
{
    int field;
    for (field = 0; field < 12; ++field) {
        if (sam_buffer_append_char(buffer, ',') < 0
                || sam_buffer_append_unsigned(
                    buffer,
                    load_u16_le_bytes(summary + field * 2)
                ) < 0) {
            return -1;
        }
    }
    return 0;
}


static int
append_columnar_annotation_tags(SamBuffer *buffer, const MateAlignment *alignment)
{
    Py_ssize_t offset;
    if (alignment->base_state_codes == NULL || alignment->read_summary == NULL) {
        PyErr_SetString(PyExc_ValueError, "SAM record has no ZT/ZR tags");
        return -1;
    }
    if (sam_buffer_append_literal(buffer, "\tzt:Z:") < 0) {
        return -1;
    }
    if (alignment->reverse_complement) {
        for (offset = alignment->query_length; offset > 0; --offset) {
            if (sam_buffer_append_char(
                    buffer, alignment->base_state_codes[offset - 1]
                ) < 0) {
                return -1;
            }
        }
    } else if (sam_buffer_append(
            buffer, alignment->base_state_codes, alignment->query_length
        ) < 0) {
        return -1;
    }
    if (sam_buffer_append_literal(buffer, "\tzr:B:S") < 0) {
        return -1;
    }
    if (append_u16_summary_values(buffer, alignment->read_summary) < 0) {
        return -1;
    }
    if (alignment->fragment_summary != NULL) {
        if (sam_buffer_append_literal(buffer, "\tzf:B:S") < 0
                || append_u16_summary_values(
                    buffer, alignment->fragment_summary
                ) < 0) {
            return -1;
        }
    }
    if (alignment->fragment_realization != NULL) {
        if (sam_buffer_append_literal(buffer, "\tzx:Z:") < 0
                || sam_buffer_append(
                    buffer,
                    alignment->fragment_realization,
                    alignment->fragment_realization_length
                ) < 0) {
            return -1;
        }
    }
    return 0;
}


static int
append_record(
    SamBuffer *buffer,
    PyObject *contig,
    PyObject *read_group_id,
    long long fragment_start,
    long long fragment_end,
    long long fragment_ordinal,
    const MateAlignment *alignment,
    const MateAlignment *other,
    long long template_length,
    int paired_end,
    PyObject *annotation_tag_suffix
)
{
    long long flag = 0;
    if (paired_end) {
        flag = 0x1 | 0x2 | (alignment->mate_index == 0 ? 0x40 : 0x80);
        if (alignment->reverse_complement) {
            flag |= 0x10;
        }
        if (other->reverse_complement) {
            flag |= 0x20;
        }
    } else if (alignment->reverse_complement) {
        flag = 0x10;
    }

    if (append_query_name(
            buffer,
            contig,
            fragment_start,
            fragment_end,
            fragment_ordinal
        ) < 0
            || sam_buffer_append_char(buffer, '\t') < 0
            || sam_buffer_append_signed(buffer, flag) < 0
            || sam_buffer_append_char(buffer, '\t') < 0
            || sam_buffer_append_unicode(buffer, contig, "contig name") < 0
            || sam_buffer_append_char(buffer, '\t') < 0
            || sam_buffer_append_signed(buffer, alignment->position) < 0
            || sam_buffer_append_literal(buffer, "\t60\t") < 0
            || sam_buffer_append(
                buffer, alignment->cigar.data, alignment->cigar.length
            ) < 0
            || sam_buffer_append_char(buffer, '\t') < 0) {
        return -1;
    }
    if (paired_end) {
        if (sam_buffer_append_literal(buffer, "=\t") < 0
                || sam_buffer_append_signed(buffer, other->position) < 0) {
            return -1;
        }
    } else if (sam_buffer_append_literal(buffer, "*\t0") < 0) {
        return -1;
    }
    if (sam_buffer_append_char(buffer, '\t') < 0
            || sam_buffer_append_signed(buffer, template_length) < 0
            || sam_buffer_append_char(buffer, '\t') < 0
            || append_reference_forward_sequence(buffer, alignment) < 0
            || sam_buffer_append_char(buffer, '\t') < 0
            || append_reference_forward_quality(buffer, alignment) < 0
            || sam_buffer_append_literal(buffer, "\tRG:Z:") < 0
            || sam_buffer_append_unicode(buffer, read_group_id, "read group") < 0
            || sam_buffer_append_literal(buffer, "\tAS:i:") < 0
            || sam_buffer_append_unsigned(
                buffer, (unsigned long long)alignment->query_length
            ) < 0) {
        return -1;
    }
    if (paired_end && sam_buffer_append_literal(buffer, "\tMQ:i:60") < 0) {
        return -1;
    }
    if (paired_end && (
            sam_buffer_append_literal(buffer, "\tMC:Z:") < 0
            || sam_buffer_append(buffer, other->cigar.data, other->cigar.length) < 0
        )) {
        return -1;
    }
    if (annotation_tag_suffix != NULL) {
        if (!PyBytes_Check(annotation_tag_suffix)) {
            PyErr_SetString(PyExc_TypeError, "details SAM tag suffix must be bytes");
            return -1;
        }
        if (sam_buffer_append(
                buffer,
                PyBytes_AS_STRING(annotation_tag_suffix),
                PyBytes_GET_SIZE(annotation_tag_suffix)
            ) < 0) {
            return -1;
        }
    } else if (append_columnar_annotation_tags(buffer, alignment) < 0) {
        return -1;
    }
    return sam_buffer_append_char(buffer, '\n');
}


static int
append_fragment_records(
    SamBuffer *buffer,
    PyObject *record_lengths,
    Py_ssize_t *record_index,
    PyObject *fragment,
    PyObject *read_group_id,
    long long contig_length,
    int paired_end,
    PyObject *annotation_tag_values
)
{
    PyObject *contig = NULL;
    PyObject *mate_values = NULL;
    PyObject *mates = NULL;
    PyObject *ordered_mates[2] = {NULL, NULL};
    MateAlignment alignments[2];
    Py_ssize_t expected_mates = paired_end ? 2 : 1;
    Py_ssize_t index;
    long long fragment_start;
    long long fragment_end;
    long long fragment_ordinal;
    long long template_lengths[2] = {0, 0};
    int result = -1;

    memset(alignments, 0, sizeof(alignments));
    contig = PyObject_GetAttrString(fragment, "contig_name");
    if (contig == NULL
            || get_integer_attribute(
                fragment, "reference_start", &fragment_start
            ) < 0
            || get_integer_attribute(fragment, "reference_end", &fragment_end) < 0
            || get_integer_attribute(
                fragment, "fragment_ordinal", &fragment_ordinal
            ) < 0) {
        goto done;
    }
    if (fragment_ordinal < 0) {
        PyErr_SetString(PyExc_ValueError, "details SAM fragment ordinal is negative");
        goto done;
    }
    mate_values = PyObject_GetAttrString(fragment, "mates");
    if (mate_values == NULL) {
        goto done;
    }
    mates = PySequence_Fast(mate_values, "details SAM mates must be a sequence");
    if (mates == NULL) {
        goto done;
    }
    if (PySequence_Fast_GET_SIZE(mates) != expected_mates) {
        PyErr_SetString(PyExc_ValueError, "details SAM mate cardinality is invalid");
        goto done;
    }
    for (index = 0; index < expected_mates; ++index) {
        PyObject *mate = PySequence_Fast_GET_ITEM(mates, index);
        long long mate_index;
        if (get_integer_attribute(mate, "mate_index", &mate_index) < 0) {
            goto done;
        }
        if (mate_index < 0 || mate_index >= expected_mates
                || ordered_mates[mate_index] != NULL) {
            PyErr_SetString(PyExc_ValueError, "details SAM mate indices are invalid");
            goto done;
        }
        ordered_mates[mate_index] = mate;
    }
    for (index = 0; index < expected_mates; ++index) {
        if (ordered_mates[index] == NULL
                || build_mate_alignment(
                    ordered_mates[index], contig_length, &alignments[index]
                ) < 0) {
            goto done;
        }
    }
    if (paired_end) {
        const long long minimum_start = alignments[0].reference_start
            < alignments[1].reference_start
            ? alignments[0].reference_start
            : alignments[1].reference_start;
        const long long maximum_end = alignments[0].reference_end
            > alignments[1].reference_end
            ? alignments[0].reference_end
            : alignments[1].reference_end;
        const long long span = maximum_end - minimum_start;
        int left_index;
        if (span <= 0 || span > INT32_MAX) {
            PyErr_SetString(PyExc_ValueError, "details SAM template length is invalid");
            goto done;
        }
        left_index = alignments[0].reference_start < alignments[1].reference_start
            || (alignments[0].reference_start == alignments[1].reference_start
                && alignments[0].mate_index < alignments[1].mate_index)
            ? 0
            : 1;
        template_lengths[left_index] = span;
        template_lengths[1 - left_index] = -span;
    }
    for (index = 0; index < expected_mates; ++index) {
        const Py_ssize_t start = buffer->length;
        PyObject *length_value;
        if (append_record(
                buffer,
                contig,
                read_group_id,
                fragment_start,
                fragment_end,
                fragment_ordinal,
                &alignments[index],
                paired_end ? &alignments[1 - index] : NULL,
                template_lengths[index],
                paired_end,
                PyTuple_GET_ITEM(annotation_tag_values, *record_index)
            ) < 0) {
            goto done;
        }
        length_value = PyLong_FromSsize_t(buffer->length - start);
        if (length_value == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(record_lengths, *record_index, length_value);
        *record_index += 1;
    }
    result = 0;

done:
    for (index = 0; index < 2; ++index) {
        mate_alignment_clear(&alignments[index]);
    }
    Py_XDECREF(mates);
    Py_XDECREF(mate_values);
    Py_XDECREF(contig);
    return result;
}


PyObject *
bsreadsim_native_format_sam_batch(PyObject *self, PyObject *args)
{
    PyObject *fragment_values;
    PyObject *read_group_id;
    PyObject *contig_length_values;
    PyObject *annotation_tag_values;
    PyObject *fragments = NULL;
    PyObject *contig_lengths = NULL;
    PyObject *record_lengths = NULL;
    PyObject *payload = NULL;
    PyObject *result = NULL;
    Py_ssize_t fragment_count;
    Py_ssize_t expected_records;
    Py_ssize_t fragment_index;
    Py_ssize_t record_index = 0;
    int paired_end;
    SamBuffer buffer = {NULL, 0, 0};

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "OpOOO:format_sam_batch",
            &fragment_values,
            &paired_end,
            &read_group_id,
            &contig_length_values,
            &annotation_tag_values
        )) {
        return NULL;
    }
    if (!PyUnicode_Check(read_group_id)) {
        PyErr_SetString(PyExc_TypeError, "details SAM read group must be str");
        return NULL;
    }
    fragments = PySequence_Fast(
        fragment_values, "details SAM fragments must be a sequence"
    );
    contig_lengths = PySequence_Fast(
        contig_length_values, "details SAM contig lengths must be a sequence"
    );
    if (fragments == NULL || contig_lengths == NULL) {
        goto done;
    }
    fragment_count = PySequence_Fast_GET_SIZE(fragments);
    if (PySequence_Fast_GET_SIZE(contig_lengths) != fragment_count) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM contig lengths disagree with fragment count"
        );
        goto done;
    }
    if (fragment_count > PY_SSIZE_T_MAX / (paired_end ? 2 : 1)) {
        PyErr_NoMemory();
        goto done;
    }
    expected_records = fragment_count * (paired_end ? 2 : 1);
    if (!PyTuple_Check(annotation_tag_values)
            || PyTuple_GET_SIZE(annotation_tag_values) != expected_records) {
        PyErr_SetString(PyExc_ValueError, "details SAM tag count is invalid");
        goto done;
    }
    record_lengths = PyTuple_New(expected_records);
    if (record_lengths == NULL) {
        goto done;
    }
    for (fragment_index = 0; fragment_index < fragment_count; ++fragment_index) {
        PyObject *length_value = PySequence_Fast_GET_ITEM(
            contig_lengths, fragment_index
        );
        long long contig_length;
        if (!PyLong_Check(length_value)
                || length_value == Py_True
                || length_value == Py_False) {
            PyErr_SetString(PyExc_TypeError, "details SAM contig length must be int");
            goto done;
        }
        contig_length = PyLong_AsLongLong(length_value);
        if ((contig_length == -1 && PyErr_Occurred()) || contig_length <= 0) {
            if (!PyErr_Occurred()) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "details SAM contig length must be positive"
                );
            }
            goto done;
        }
        if (append_fragment_records(
                &buffer,
                record_lengths,
                &record_index,
                PySequence_Fast_GET_ITEM(fragments, fragment_index),
                read_group_id,
                contig_length,
                paired_end,
                annotation_tag_values
            ) < 0) {
            goto done;
        }
    }
    payload = PyBytes_FromStringAndSize(buffer.data, buffer.length);
    if (payload == NULL) {
        goto done;
    }
    result = PyTuple_New(2);
    if (result == NULL) {
        goto done;
    }
    PyTuple_SET_ITEM(result, 0, payload);
    payload = NULL;
    PyTuple_SET_ITEM(result, 1, record_lengths);
    record_lengths = NULL;

done:
    sam_buffer_clear(&buffer);
    Py_XDECREF(payload);
    Py_XDECREF(record_lengths);
    Py_XDECREF(contig_lengths);
    Py_XDECREF(fragments);
    return result;
}


#define SAM_COLUMN_COUNT 15


typedef struct {
    Py_buffer view;
    int acquired;
} SamColumn;


static uint32_t
load_u32_le(const SamColumn *column, Py_ssize_t index)
{
    const unsigned char *value =
        (const unsigned char *)column->view.buf + index * 4;
    return (uint32_t)value[0]
        | ((uint32_t)value[1] << 8)
        | ((uint32_t)value[2] << 16)
        | ((uint32_t)value[3] << 24);
}


static void
sam_columns_clear(SamColumn *columns)
{
    int index;
    for (index = 0; index < SAM_COLUMN_COUNT; ++index) {
        if (columns[index].acquired) {
            PyBuffer_Release(&columns[index].view);
            columns[index].acquired = 0;
        }
    }
}


static int
sam_columns_acquire(PyObject *values, SamColumn *columns)
{
    Py_ssize_t index;
    if (!PyTuple_Check(values)
            || PyTuple_GET_SIZE(values) != SAM_COLUMN_COUNT) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM column tuple has the wrong cardinality"
        );
        return -1;
    }
    for (index = 0; index < SAM_COLUMN_COUNT; ++index) {
        if (PyObject_GetBuffer(
                PyTuple_GET_ITEM(values, index),
                &columns[index].view,
                PyBUF_SIMPLE
            ) < 0) {
            sam_columns_clear(columns);
            return -1;
        }
        columns[index].acquired = 1;
    }
    return 0;
}


static int
require_column_size(
    const SamColumn *column,
    Py_ssize_t count,
    Py_ssize_t width,
    const char *name
)
{
    if (count < 0 || width <= 0 || count > PY_SSIZE_T_MAX / width
            || column->view.len != count * width) {
        PyErr_Format(PyExc_ValueError, "details SAM %s has the wrong size", name);
        return -1;
    }
    return 0;
}


static int
build_columnar_alignment(
    const SamColumn *columns,
    Py_ssize_t fragment_index,
    Py_ssize_t mate_row,
    long long contig_length,
    const unsigned char *sequences,
    const unsigned char *base_state_codes,
    const unsigned char *read_summaries,
    Py_ssize_t read_length,
    int quality_byte,
    MateAlignment *alignment
)
{
    const uint32_t mate_template_begin = load_u32_le(&columns[3], mate_row);
    const uint32_t mate_template_end = load_u32_le(&columns[4], mate_row);
    const uint32_t projection_begin = load_u32_le(&columns[7], fragment_index);
    const uint32_t projection_end = load_u32_le(&columns[7], fragment_index + 1);
    const uint32_t event_begin = load_u32_le(&columns[11], fragment_index);
    const uint32_t event_end = load_u32_le(&columns[11], fragment_index + 1);
    uint32_t projection_index;
    uint32_t cursor = mate_template_begin;
    int have_mapped = 0;
    long long previous_reference_end = 0;
    long long mapped_start = 0;
    long long mapped_end = 0;
    uint64_t query_consumed = 0;
    char pending_operation = '\0';
    uint64_t pending_length = 0;

    memset(alignment, 0, sizeof(*alignment));
    if (mate_template_end <= mate_template_begin
            || (Py_ssize_t)(mate_template_end - mate_template_begin) != read_length) {
        PyErr_SetString(PyExc_ValueError, "details SAM mate template slice is invalid");
        return -1;
    }
    alignment->mate_index = ((const unsigned char *)columns[5].view.buf)[mate_row];
    alignment->reverse_complement =
        ((const unsigned char *)columns[6].view.buf)[mate_row] != 0;
    alignment->sequence = (const char *)sequences + mate_row * read_length;
    alignment->base_state_codes = (const char *)base_state_codes + mate_row * read_length;
    alignment->read_summary = read_summaries + mate_row * 12 * 2;
    alignment->quality = NULL;
    alignment->query_length = read_length;
    alignment->quality_byte = quality_byte;

    for (projection_index = projection_begin;
            projection_index < projection_end;
            ++projection_index) {
        const uint32_t segment_template_begin =
            load_u32_le(&columns[8], projection_index);
        const uint32_t segment_template_end =
            load_u32_le(&columns[9], projection_index);
        const uint32_t segment_reference_begin =
            load_u32_le(&columns[10], projection_index);
        const uint32_t intersection_begin = segment_template_begin > mate_template_begin
            ? segment_template_begin
            : mate_template_begin;
        const uint32_t intersection_end = segment_template_end < mate_template_end
            ? segment_template_end
            : mate_template_end;
        long long reference_begin;
        uint64_t mapped_length;
        if (intersection_begin >= intersection_end) {
            continue;
        }
        if (intersection_begin > cursor) {
            const uint64_t insertion_length = intersection_begin - cursor;
            if (push_cigar_operation(
                    &alignment->cigar,
                    &pending_operation,
                    &pending_length,
                    'I',
                    insertion_length
                ) < 0) {
                goto error;
            }
            query_consumed += insertion_length;
        }
        reference_begin = (long long)segment_reference_begin
            + (intersection_begin - segment_template_begin);
        if (reference_begin < 0 || reference_begin >= contig_length) {
            PyErr_SetString(
                PyExc_ValueError,
                "details SAM projection is outside its contig"
            );
            goto error;
        }
        if (have_mapped) {
            const long long deletion_length =
                reference_begin - previous_reference_end;
            if (deletion_length < 0) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "details SAM projection references are not increasing"
                );
                goto error;
            }
            if (deletion_length > 0 && push_cigar_operation(
                    &alignment->cigar,
                    &pending_operation,
                    &pending_length,
                    'D',
                    (uint64_t)deletion_length
                ) < 0) {
                goto error;
            }
        } else {
            mapped_start = reference_begin;
            have_mapped = 1;
        }
        mapped_length = intersection_end - intersection_begin;
        if (push_cigar_operation(
                &alignment->cigar,
                &pending_operation,
                &pending_length,
                'M',
                mapped_length
            ) < 0) {
            goto error;
        }
        query_consumed += mapped_length;
        mapped_end = reference_begin + (long long)mapped_length;
        previous_reference_end = mapped_end;
        cursor = intersection_end;
    }
    if (cursor < mate_template_end) {
        const uint64_t insertion_length = mate_template_end - cursor;
        if (push_cigar_operation(
                &alignment->cigar,
                &pending_operation,
                &pending_length,
                'I',
                insertion_length
            ) < 0) {
            goto error;
        }
        query_consumed += insertion_length;
    }
    if (flush_cigar_operation(
            &alignment->cigar, &pending_operation, &pending_length
        ) < 0) {
        goto error;
    }
    if (query_consumed != (uint64_t)read_length) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM columnar CIGAR does not consume the complete read"
        );
        goto error;
    }
    if (have_mapped) {
        alignment->reference_start = mapped_start;
        alignment->reference_end = mapped_end;
    } else {
        uint32_t variant_index;
        int have_anchor = 0;
        long long anchor = 0;
        for (variant_index = event_begin; variant_index < event_end; ++variant_index) {
            const uint32_t template_begin = load_u32_le(&columns[13], variant_index);
            const uint32_t template_end = load_u32_le(&columns[14], variant_index);
            const long long candidate = load_u32_le(&columns[12], variant_index);
            if (template_begin >= mate_template_end
                    || template_end <= mate_template_begin) {
                continue;
            }
            if (have_anchor && candidate != anchor) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "details SAM insertion-only mate has ambiguous anchors"
                );
                goto error;
            }
            anchor = candidate;
            have_anchor = 1;
        }
        if (!have_anchor || anchor < 0 || anchor > contig_length) {
            PyErr_SetString(
                PyExc_ValueError,
                "details SAM insertion-only mate has no valid anchor"
            );
            goto error;
        }
        alignment->reference_start = anchor < contig_length
            ? anchor
            : contig_length - 1;
        alignment->reference_end = alignment->reference_start + 1;
    }
    alignment->position = alignment->reference_start + 1;
    return 0;

error:
    mate_alignment_clear(alignment);
    return -1;
}


PyObject *
bsreadsim_native_format_sam_columns(PyObject *self, PyObject *args)
{
    PyObject *contig_names;
    PyObject *column_values;
    PyObject *read_group_id;
    PyObject *contig_length_values;
    PyObject *contig_lengths = NULL;
    Py_buffer sequences = {0};
    Py_buffer base_state_codes = {0};
    Py_buffer read_summaries = {0};
    Py_buffer fragment_summaries = {0};
    PyObject *fragment_summaries_value;
    PyObject *fragment_realizations_value;
    SamColumn columns[SAM_COLUMN_COUNT] = {0};
    PyObject *record_lengths = NULL;
    PyObject *payload = NULL;
    PyObject *result = NULL;
    SamBuffer buffer = {NULL, 0, 0};
    Py_ssize_t fragment_count;
    Py_ssize_t mate_count;
    Py_ssize_t expected_records;
    Py_ssize_t fragment_index;
    Py_ssize_t record_index = 0;
    Py_ssize_t read_length;
    unsigned long long first_fragment_ordinal;
    int quality_byte;
    int paired_end;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "OOy*y*y*OOniiKOO:format_sam_columns",
            &contig_names,
            &column_values,
            &sequences,
            &base_state_codes,
            &read_summaries,
            &fragment_summaries_value,
            &fragment_realizations_value,
            &read_length,
            &quality_byte,
            &paired_end,
            &first_fragment_ordinal,
            &read_group_id,
            &contig_length_values
        )) {
        return NULL;
    }
    if (fragment_summaries_value != Py_None
            && PyObject_GetBuffer(
                fragment_summaries_value,
                &fragment_summaries,
                PyBUF_SIMPLE
            ) < 0) {
        goto done;
    }
    if (!PyTuple_Check(contig_names)) {
        PyErr_SetString(PyExc_TypeError, "details SAM contig names must be a tuple");
        goto done;
    }
    if (!PyUnicode_Check(read_group_id)) {
        PyErr_SetString(PyExc_TypeError, "details SAM read group must be str");
        goto done;
    }
    if (read_length <= 0 || quality_byte < 33 || quality_byte > 126) {
        PyErr_SetString(PyExc_ValueError, "details SAM read geometry is invalid");
        goto done;
    }
    if (sam_columns_acquire(column_values, columns) < 0) {
        goto done;
    }
    fragment_count = PyTuple_GET_SIZE(contig_names);
    if (fragment_realizations_value != Py_None
            && (!PyTuple_Check(fragment_realizations_value)
                || PyTuple_GET_SIZE(fragment_realizations_value)
                    != fragment_count)) {
        PyErr_SetString(
            PyExc_ValueError,
            "details SAM fragment realizations disagree with fragments"
        );
        goto done;
    }
    contig_lengths = PySequence_Fast(
        contig_length_values, "details SAM contig lengths must be a sequence"
    );
    if (contig_lengths == NULL
            || PySequence_Fast_GET_SIZE(contig_lengths) != fragment_count) {
        if (!PyErr_Occurred()) {
            PyErr_SetString(
                PyExc_ValueError,
                "details SAM contig lengths disagree with fragments"
            );
        }
        goto done;
    }
    if (require_column_size(&columns[0], fragment_count, 4, "reference begins") < 0
            || require_column_size(
                &columns[1], fragment_count, 4, "reference ends"
            ) < 0
            || require_column_size(
                &columns[2], fragment_count + 1, 4, "mate offsets"
            ) < 0
            || require_column_size(
                &columns[7], fragment_count + 1, 4, "projection offsets"
            ) < 0
            || require_column_size(
                &columns[11], fragment_count + 1, 4, "event offsets"
            ) < 0) {
        goto done;
    }
    mate_count = load_u32_le(&columns[2], fragment_count);
    expected_records = fragment_count * (paired_end ? 2 : 1);
    if (mate_count != expected_records
            || require_column_size(&columns[3], mate_count, 4, "mate starts") < 0
            || require_column_size(&columns[4], mate_count, 4, "mate ends") < 0
            || require_column_size(&columns[5], mate_count, 1, "mate indices") < 0
            || require_column_size(&columns[6], mate_count, 1, "mate strands") < 0
            || sequences.len != mate_count * read_length
            || base_state_codes.len != mate_count * read_length
            || read_summaries.len != mate_count * 12 * 2
            || (fragment_summaries.obj != NULL
                && fragment_summaries.len != fragment_count * 12 * 2)) {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError, "details SAM mate columns are invalid");
        }
        goto done;
    }
    {
        const Py_ssize_t projection_count =
            load_u32_le(&columns[7], fragment_count);
        const Py_ssize_t variant_count = load_u32_le(&columns[11], fragment_count);
        if (require_column_size(
                &columns[8], projection_count, 4, "projection starts"
            ) < 0
                || require_column_size(
                    &columns[9], projection_count, 4, "projection ends"
                ) < 0
                || require_column_size(
                    &columns[10], projection_count, 4, "projection references"
                ) < 0
                || require_column_size(
                    &columns[12], variant_count, 4, "event references"
                ) < 0
                || require_column_size(
                    &columns[13], variant_count, 4, "event starts"
                ) < 0
                || require_column_size(
                    &columns[14], variant_count, 4, "event ends"
                ) < 0) {
            goto done;
        }
    }
    record_lengths = PyTuple_New(expected_records);
    if (record_lengths == NULL) {
        goto done;
    }
    for (fragment_index = 0; fragment_index < fragment_count; ++fragment_index) {
        const uint32_t mate_begin = load_u32_le(&columns[2], fragment_index);
        const uint32_t mate_end = load_u32_le(&columns[2], fragment_index + 1);
        const Py_ssize_t expected_mates = paired_end ? 2 : 1;
        PyObject *contig = PyTuple_GET_ITEM(contig_names, fragment_index);
        PyObject *length_value = PySequence_Fast_GET_ITEM(
            contig_lengths, fragment_index
        );
        MateAlignment alignments[2];
        Py_ssize_t index;
        long long contig_length;
        long long template_lengths[2] = {0, 0};
        const char *fragment_realization = NULL;
        Py_ssize_t fragment_realization_length = 0;
        Py_ssize_t realization_offset;
        memset(alignments, 0, sizeof(alignments));
        if (fragment_realizations_value != Py_None) {
            PyObject *value = PyTuple_GET_ITEM(
                fragment_realizations_value, fragment_index
            );
            if (!PyBytes_Check(value)) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "details SAM fragment realization must be bytes"
                );
                goto fragment_error;
            }
            fragment_realization = PyBytes_AS_STRING(value);
            fragment_realization_length = PyBytes_GET_SIZE(value);
            if (fragment_realization_length <= 0) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "details SAM fragment realization must not be empty"
                );
                goto fragment_error;
            }
            for (realization_offset = 0;
                    realization_offset < fragment_realization_length;
                    ++realization_offset) {
                const unsigned char byte = (unsigned char)
                    fragment_realization[realization_offset];
                if (byte < 32 || byte > 126 || byte == '\t') {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "details SAM fragment realization must be printable ASCII"
                    );
                    goto fragment_error;
                }
            }
        }
        if (!PyLong_Check(length_value)
                || length_value == Py_True
                || length_value == Py_False
                || (contig_length = PyLong_AsLongLong(length_value)) <= 0) {
            if (!PyErr_Occurred()) {
                PyErr_SetString(PyExc_ValueError, "details SAM contig length is invalid");
            }
            goto fragment_error;
        }
        if (mate_end - mate_begin != expected_mates) {
            PyErr_SetString(PyExc_ValueError, "details SAM mate count is invalid");
            goto fragment_error;
        }
        for (index = 0; index < expected_mates; ++index) {
            const Py_ssize_t mate_row = mate_begin + index;
            const unsigned char mate_index =
                ((const unsigned char *)columns[5].view.buf)[mate_row];
            if (mate_index >= expected_mates
                    || build_columnar_alignment(
                        columns,
                        fragment_index,
                        mate_row,
                        contig_length,
                        (const unsigned char *)sequences.buf,
                        (const unsigned char *)base_state_codes.buf,
                        (const unsigned char *)read_summaries.buf,
                        read_length,
                        quality_byte,
                        &alignments[mate_index]
                    ) < 0) {
                goto fragment_error;
            }
            alignments[mate_index].fragment_summary =
                fragment_summaries.obj == NULL
                ? NULL
                : (const unsigned char *)fragment_summaries.buf
                    + fragment_index * 12 * 2;
            alignments[mate_index].fragment_realization = fragment_realization;
            alignments[mate_index].fragment_realization_length =
                fragment_realization_length;
        }
        if (paired_end) {
            const long long minimum_start = alignments[0].reference_start
                < alignments[1].reference_start
                ? alignments[0].reference_start
                : alignments[1].reference_start;
            const long long maximum_end = alignments[0].reference_end
                > alignments[1].reference_end
                ? alignments[0].reference_end
                : alignments[1].reference_end;
            const long long span = maximum_end - minimum_start;
            const int left_index =
                alignments[0].reference_start < alignments[1].reference_start
                || (alignments[0].reference_start == alignments[1].reference_start
                    && alignments[0].mate_index < alignments[1].mate_index)
                ? 0
                : 1;
            if (span <= 0 || span > INT32_MAX) {
                PyErr_SetString(PyExc_ValueError, "details SAM template length is invalid");
                goto fragment_error;
            }
            template_lengths[left_index] = span;
            template_lengths[1 - left_index] = -span;
        }
        for (index = 0; index < expected_mates; ++index) {
            const Py_ssize_t start = buffer.length;
            PyObject *record_length;
            const unsigned long long ordinal =
                first_fragment_ordinal + (unsigned long long)fragment_index;
            if (ordinal > (unsigned long long)INT64_MAX
                    || append_record(
                        &buffer,
                        contig,
                        read_group_id,
                        load_u32_le(&columns[0], fragment_index),
                        load_u32_le(&columns[1], fragment_index),
                        (long long)ordinal,
                        &alignments[index],
                        paired_end ? &alignments[1 - index] : NULL,
                        template_lengths[index],
                        paired_end,
                        NULL
                    ) < 0) {
                if (!PyErr_Occurred()) {
                    PyErr_SetString(
                        PyExc_OverflowError,
                        "details SAM fragment ordinal exceeds int64"
                    );
                }
                goto fragment_error;
            }
            record_length = PyLong_FromSsize_t(buffer.length - start);
            if (record_length == NULL) {
                goto fragment_error;
            }
            PyTuple_SET_ITEM(record_lengths, record_index, record_length);
            record_index += 1;
        }
        for (index = 0; index < 2; ++index) {
            mate_alignment_clear(&alignments[index]);
        }
        continue;

fragment_error:
        for (index = 0; index < 2; ++index) {
            mate_alignment_clear(&alignments[index]);
        }
        goto done;
    }
    payload = PyBytes_FromStringAndSize(buffer.data, buffer.length);
    if (payload == NULL) {
        goto done;
    }
    result = PyTuple_New(2);
    if (result == NULL) {
        goto done;
    }
    PyTuple_SET_ITEM(result, 0, payload);
    payload = NULL;
    PyTuple_SET_ITEM(result, 1, record_lengths);
    record_lengths = NULL;

done:
    sam_buffer_clear(&buffer);
    sam_columns_clear(columns);
    if (sequences.obj != NULL) {
        PyBuffer_Release(&sequences);
    }
    if (base_state_codes.obj != NULL) {
        PyBuffer_Release(&base_state_codes);
    }
    if (read_summaries.obj != NULL) {
        PyBuffer_Release(&read_summaries);
    }
    if (fragment_summaries.obj != NULL) {
        PyBuffer_Release(&fragment_summaries);
    }
    Py_XDECREF(contig_lengths);
    Py_XDECREF(payload);
    Py_XDECREF(record_lengths);
    return result;
}
