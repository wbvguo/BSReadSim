#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "api.h"


#define COMMON_COLUMN_COUNT 18
#define DETAIL_COLUMN_COUNT 22
#define NO_REFERENCE_POSITION UINT32_C(0xffffffff)


typedef struct {
    Py_buffer view;
    int acquired;
} ColumnBuffer;


static const char *common_column_names[COMMON_COLUMN_COUNT] = {
    "batch.contig_indices",
    "batch.reference_starts",
    "batch.reference_ends",
    "batch.template_offsets",
    "batch.mate_offsets",
    "batch.site_offsets",
    "batch.mate_template_starts",
    "batch.mate_template_ends",
    "batch.site_template_offsets",
    "batch.site_probabilities",
    "batch.haplotypes",
    "batch.capture_strands",
    "batch.mate_indices",
    "batch.mate_reverse_complements",
    "batch.site_contexts",
    "batch.methylation_sources",
    "batch.site_alleles",
    "batch.template_bases"
};


static const char *detail_column_names[DETAIL_COLUMN_COUNT] = {
    "details.projection_offsets",
    "details.variant_offsets",
    "details.original_n_offsets",
    "details.projection_template_starts",
    "details.projection_template_ends",
    "details.projection_reference_starts",
    "details.variant_indices",
    "details.variant_id_offsets",
    "details.variant_reference_starts",
    "details.variant_reference_ends",
    "details.variant_template_starts",
    "details.variant_template_ends",
    "details.variant_ref_offsets",
    "details.variant_alt_offsets",
    "details.site_reference_positions",
    "details.original_n_template_offsets",
    "details.variant_sources",
    "details.variant_kinds",
    "details.variant_phased_haplotypes",
    "details.variant_ids",
    "details.variant_ref_bases",
    "details.variant_alt_bases"
};


static uint32_t
load_u32_le(const ColumnBuffer *column, uint32_t index)
{
    const unsigned char *data = (const unsigned char *)column->view.buf;
    const size_t offset = (size_t)index * 4U;
    return
        (uint32_t)data[offset]
        | ((uint32_t)data[offset + 1U] << 8)
        | ((uint32_t)data[offset + 2U] << 16)
        | ((uint32_t)data[offset + 3U] << 24);
}


static float
load_f32_le(const ColumnBuffer *column, uint32_t index)
{
    const uint32_t bits = load_u32_le(column, index);
    float value;
    memcpy(&value, &bits, sizeof(value));
    return value;
}


static unsigned char
load_u8(const ColumnBuffer *column, uint32_t index)
{
    return ((const unsigned char *)column->view.buf)[index];
}


static void
release_columns(ColumnBuffer *columns, Py_ssize_t count)
{
    Py_ssize_t index;
    for (index = 0; index < count; ++index) {
        if (columns[index].acquired) {
            PyBuffer_Release(&columns[index].view);
            columns[index].acquired = 0;
        }
    }
}


static int
acquire_columns(
    PyObject *value,
    Py_ssize_t expected_count,
    const char *const *names,
    ColumnBuffer *columns,
    const char *group_name
)
{
    Py_ssize_t index;
    if (!PyTuple_Check(value) || PyTuple_GET_SIZE(value) != expected_count) {
        PyErr_Format(
            PyExc_TypeError,
            "%s must be a %zd-column tuple",
            group_name,
            expected_count
        );
        return 0;
    }
    for (index = 0; index < expected_count; ++index) {
        PyObject *item = PyTuple_GET_ITEM(value, index);
        if (PyObject_GetBuffer(item, &columns[index].view, PyBUF_CONTIG_RO) < 0) {
            return 0;
        }
        columns[index].acquired = 1;
        if (!columns[index].view.readonly) {
            PyErr_Format(
                PyExc_BufferError,
                "%s must expose an immutable buffer",
                names[index]
            );
            return 0;
        }
        if (columns[index].view.len < 0) {
            PyErr_Format(PyExc_BufferError, "%s has an invalid size", names[index]);
            return 0;
        }
    }
    return 1;
}


static int
parse_u32(PyObject *value, const char *name, uint32_t *result)
{
    unsigned long long converted;
    if (PyBool_Check(value) || !PyLong_Check(value)) {
        PyErr_Format(PyExc_TypeError, "%s must be an integer", name);
        return 0;
    }
    converted = PyLong_AsUnsignedLongLong(value);
    if (converted == (unsigned long long)-1 && PyErr_Occurred()) {
        return 0;
    }
    if (converted > UINT32_MAX) {
        PyErr_Format(PyExc_OverflowError, "%s exceeds u32", name);
        return 0;
    }
    *result = (uint32_t)converted;
    return 1;
}


static int
derive_u32_count(
    const ColumnBuffer *column,
    uint32_t item_size,
    const char *name,
    uint32_t *result
)
{
    const uint64_t length = (uint64_t)column->view.len;
    const uint64_t count = length / item_size;
    if (length % item_size != 0) {
        PyErr_Format(PyExc_ValueError, "%s byte length is not canonical", name);
        return 0;
    }
    if (count > UINT32_MAX) {
        PyErr_Format(PyExc_OverflowError, "%s row count exceeds u32", name);
        return 0;
    }
    *result = (uint32_t)count;
    return 1;
}


static int
require_column_size(
    const ColumnBuffer *column,
    uint64_t count,
    uint32_t item_size,
    const char *name
)
{
    const uint64_t expected = count * item_size;
    if ((uint64_t)column->view.len != expected) {
        PyErr_Format(
            PyExc_ValueError,
            "%s has %zd bytes; expected %llu",
            name,
            column->view.len,
            (unsigned long long)expected
        );
        return 0;
    }
    return 1;
}


static int
validate_prefix(
    const ColumnBuffer *column,
    uint32_t row_count,
    uint32_t flat_count,
    const char *name
)
{
    uint32_t index;
    uint32_t previous;
    if (!require_column_size(column, (uint64_t)row_count + 1U, 4, name)) {
        return 0;
    }
    previous = load_u32_le(column, 0);
    if (previous != 0) {
        PyErr_Format(PyExc_ValueError, "%s must begin at zero", name);
        return 0;
    }
    for (index = 1; index <= row_count; ++index) {
        const uint32_t current = load_u32_le(column, index);
        if (current < previous) {
            PyErr_Format(PyExc_ValueError, "%s must be monotone", name);
            return 0;
        }
        previous = current;
    }
    if (previous != flat_count) {
        PyErr_Format(
            PyExc_ValueError,
            "%s must end at %u",
            name,
            (unsigned int)flat_count
        );
        return 0;
    }
    return 1;
}


static int
is_context(unsigned char value)
{
    return
        value == 1 || value == 3 || value == 7
        || value == 9 || value == 11 || value == 15;
}


static int
is_cytosine_context(unsigned char value)
{
    return value == 1 || value == 3 || value == 7;
}


static int
is_guanine_context(unsigned char value)
{
    return value == 9 || value == 11 || value == 15;
}


static int
validate_common_columns(
    ColumnBuffer *common,
    const uint32_t *contig_lengths,
    uint32_t contig_count,
    uint32_t first_ordinal,
    uint32_t mates_per_fragment,
    uint32_t read_length_r1,
    uint32_t read_length_r2,
    int has_expected_ordinal,
    uint32_t expected_ordinal,
    uint32_t *fragment_count_out,
    uint32_t *template_count_out,
    uint32_t *mate_count_out,
    uint32_t *site_count_out,
    uint32_t *maximum_template_length_out
)
{
    uint32_t fragment_count;
    uint32_t template_count;
    uint32_t mate_count;
    uint32_t site_count;
    uint32_t maximum_template_length = 0;
    uint32_t row;
    uint32_t index;

    if (
        !derive_u32_count(&common[0], 4, common_column_names[0], &fragment_count)
        || !derive_u32_count(&common[17], 1, common_column_names[17], &template_count)
        || !derive_u32_count(&common[12], 1, common_column_names[12], &mate_count)
        || !derive_u32_count(&common[8], 4, common_column_names[8], &site_count)
    ) {
        return 0;
    }
    if (fragment_count == 0) {
        PyErr_SetString(
            PyExc_ValueError,
            "fragment batch must contain at least one row"
        );
        return 0;
    }
    if ((uint64_t)first_ordinal + fragment_count > UINT32_MAX) {
        PyErr_SetString(
            PyExc_ValueError,
            "batch fragment ordinal range exceeds u32"
        );
        return 0;
    }
    if (has_expected_ordinal && first_ordinal != expected_ordinal) {
        PyErr_SetString(
            PyExc_ValueError,
            "fragment batch ordinal range is not consecutive"
        );
        return 0;
    }

    if (
        !require_column_size(&common[1], fragment_count, 4, common_column_names[1])
        || !require_column_size(&common[2], fragment_count, 4, common_column_names[2])
        || !validate_prefix(
            &common[3], fragment_count, template_count, common_column_names[3]
        )
        || !validate_prefix(
            &common[4], fragment_count, mate_count, common_column_names[4]
        )
        || !validate_prefix(
            &common[5], fragment_count, site_count, common_column_names[5]
        )
        || !require_column_size(&common[6], mate_count, 4, common_column_names[6])
        || !require_column_size(&common[7], mate_count, 4, common_column_names[7])
        || !require_column_size(&common[9], site_count, 4, common_column_names[9])
        || !require_column_size(&common[10], fragment_count, 1, common_column_names[10])
        || !require_column_size(&common[11], fragment_count, 1, common_column_names[11])
        || !require_column_size(&common[13], mate_count, 1, common_column_names[13])
        || !require_column_size(&common[14], site_count, 1, common_column_names[14])
        || !require_column_size(&common[15], site_count, 1, common_column_names[15])
        || !require_column_size(&common[16], site_count, 1, common_column_names[16])
    ) {
        return 0;
    }

    for (index = 0; index < template_count; ++index) {
        if (load_u8(&common[17], index) > 4) {
            PyErr_SetString(
                PyExc_ValueError,
                "batch.template_bases contains an invalid base code"
            );
            return 0;
        }
    }

    for (row = 0; row < fragment_count; ++row) {
        const uint32_t contig_index = load_u32_le(&common[0], row);
        const uint32_t reference_begin = load_u32_le(&common[1], row);
        const uint32_t reference_end = load_u32_le(&common[2], row);
        const uint32_t template_begin = load_u32_le(&common[3], row);
        const uint32_t template_end = load_u32_le(&common[3], row + 1U);
        const uint32_t template_length = template_end - template_begin;
        const uint32_t mate_begin = load_u32_le(&common[4], row);
        const uint32_t mate_end = load_u32_le(&common[4], row + 1U);
        const uint32_t site_begin = load_u32_le(&common[5], row);
        const uint32_t site_end = load_u32_le(&common[5], row + 1U);
        uint32_t local_mate;
        uint32_t previous_site = 0;
        int has_previous_site = 0;

        if (contig_index >= contig_count) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment contig index is outside the header"
            );
            return 0;
        }
        if (
            reference_begin >= reference_end
            || reference_end > contig_lengths[contig_index]
        ) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment reference envelope is invalid"
            );
            return 0;
        }
        if (load_u8(&common[10], row) > 1) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment haplotype must be zero or one"
            );
            return 0;
        }
        if (load_u8(&common[11], row) > 2) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment.capture_strand is outside CaptureStrand"
            );
            return 0;
        }
        if (template_length == 0) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment template must not be empty"
            );
            return 0;
        }
        if (template_length > maximum_template_length) {
            maximum_template_length = template_length;
        }
        if (mate_end - mate_begin != mates_per_fragment) {
            PyErr_SetString(
                PyExc_ValueError,
                "fragment mate count disagrees with header"
            );
            return 0;
        }
        for (local_mate = 0; local_mate < mates_per_fragment; ++local_mate) {
            const uint32_t mate_index = mate_begin + local_mate;
            const uint32_t begin = load_u32_le(&common[6], mate_index);
            const uint32_t end = load_u32_le(&common[7], mate_index);
            const uint32_t expected_length =
                local_mate == 0 ? read_length_r1 : read_length_r2;
            if (load_u8(&common[12], mate_index) != local_mate) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "mate rows must be ordered as 0 or 0,1"
                );
                return 0;
            }
            if (load_u8(&common[13], mate_index) > 1) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "mate reverse-complement flag must be boolean"
                );
                return 0;
            }
            if (
                begin >= end
                || end > template_length
                || end - begin != expected_length
            ) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "mate template slice disagrees with header"
                );
                return 0;
            }
        }
        for (index = site_begin; index < site_end; ++index) {
            const uint32_t local_offset = load_u32_le(&common[8], index);
            const float probability = load_f32_le(&common[9], index);
            const unsigned char context = load_u8(&common[14], index);
            const unsigned char source = load_u8(&common[15], index);
            const unsigned char allele = load_u8(&common[16], index);
            unsigned char base;
            if (local_offset >= template_length) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site offset exceeds its template"
                );
                return 0;
            }
            if (has_previous_site && local_offset <= previous_site) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site offsets must be strictly increasing"
                );
                return 0;
            }
            previous_site = local_offset;
            has_previous_site = 1;
            if (!isfinite(probability) || probability < 0.0F || probability > 1.0F) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site.probability must be finite and in [0,1]"
                );
                return 0;
            }
            if (!is_context(context)) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site.context is outside MethylationContext"
                );
                return 0;
            }
            if (source < 1 || source > 4) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site.source is outside MethylationSource"
                );
                return 0;
            }
            if (allele > 2) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site.allele is outside MethylationAllele"
                );
                return 0;
            }
            base = load_u8(&common[17], template_begin + local_offset);
            if (
                (is_cytosine_context(context) && base != 1)
                || (is_guanine_context(context) && base != 2)
            ) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site context is incompatible with template base"
                );
                return 0;
            }
        }
    }

    *fragment_count_out = fragment_count;
    *template_count_out = template_count;
    *mate_count_out = mate_count;
    *site_count_out = site_count;
    *maximum_template_length_out = maximum_template_length;
    return 1;
}


static int
validate_details(
    ColumnBuffer *common,
    ColumnBuffer *details,
    uint32_t fragment_count,
    uint32_t site_count,
    uint32_t maximum_template_length
)
{
    uint32_t projection_count;
    uint32_t variant_count;
    uint32_t variant_id_byte_count;
    uint32_t ref_base_count;
    uint32_t alt_base_count;
    uint32_t original_n_count;
    unsigned char *projection_cover = NULL;
    unsigned char *insertion_cover = NULL;
    unsigned char *event_cover = NULL;
    unsigned char *observed_n = NULL;
    uint32_t *mapped_positions = NULL;
    uint32_t row;
    uint32_t index;
    int result = 0;

    if (
        !derive_u32_count(&details[3], 4, detail_column_names[3], &projection_count)
        || !derive_u32_count(&details[6], 4, detail_column_names[6], &variant_count)
        || !derive_u32_count(&details[19], 1, detail_column_names[19], &variant_id_byte_count)
        || !derive_u32_count(&details[20], 1, detail_column_names[20], &ref_base_count)
        || !derive_u32_count(&details[21], 1, detail_column_names[21], &alt_base_count)
        || !derive_u32_count(&details[15], 4, detail_column_names[15], &original_n_count)
    ) {
        goto done;
    }
    if (
        !validate_prefix(
            &details[0], fragment_count, projection_count, detail_column_names[0]
        )
        || !validate_prefix(
            &details[1], fragment_count, variant_count, detail_column_names[1]
        )
        || !validate_prefix(
            &details[2], fragment_count, original_n_count, detail_column_names[2]
        )
        || !require_column_size(&details[4], projection_count, 4, detail_column_names[4])
        || !require_column_size(&details[5], projection_count, 4, detail_column_names[5])
        || !validate_prefix(
            &details[7], variant_count, variant_id_byte_count, detail_column_names[7]
        )
        || !require_column_size(&details[8], variant_count, 4, detail_column_names[8])
        || !require_column_size(&details[9], variant_count, 4, detail_column_names[9])
        || !require_column_size(&details[10], variant_count, 4, detail_column_names[10])
        || !require_column_size(&details[11], variant_count, 4, detail_column_names[11])
        || !validate_prefix(
            &details[12], variant_count, ref_base_count, detail_column_names[12]
        )
        || !validate_prefix(
            &details[13], variant_count, alt_base_count, detail_column_names[13]
        )
        || !require_column_size(&details[14], site_count, 4, detail_column_names[14])
        || !require_column_size(&details[16], variant_count, 1, detail_column_names[16])
        || !require_column_size(&details[17], variant_count, 1, detail_column_names[17])
        || !require_column_size(&details[18], variant_count, 1, detail_column_names[18])
    ) {
        goto done;
    }
    for (index = 0; index < variant_count; ++index) {
        const uint32_t id_begin = load_u32_le(&details[7], index);
        const uint32_t id_end = load_u32_le(&details[7], index + 1U);
        const unsigned char source = load_u8(&details[16], index);
        if (id_begin == id_end) {
            PyErr_SetString(PyExc_ValueError, "variant ID must not be empty");
            goto done;
        }
        if (source < 1 || source > 2) {
            PyErr_SetString(PyExc_ValueError, "variant source is invalid");
            goto done;
        }
    }
    for (index = 0; index < ref_base_count; ++index) {
        if (load_u8(&details[20], index) > 3) {
            PyErr_SetString(
                PyExc_ValueError,
                "details.variant_ref_bases contains an invalid base code"
            );
            goto done;
        }
    }
    for (index = 0; index < alt_base_count; ++index) {
        if (load_u8(&details[21], index) > 3) {
            PyErr_SetString(
                PyExc_ValueError,
                "details.variant_alt_bases contains an invalid base code"
            );
            goto done;
        }
    }

    projection_cover = PyMem_Malloc((size_t)maximum_template_length);
    insertion_cover = PyMem_Malloc((size_t)maximum_template_length);
    event_cover = PyMem_Malloc((size_t)maximum_template_length);
    observed_n = PyMem_Malloc((size_t)maximum_template_length);
    mapped_positions = PyMem_Malloc(
        (size_t)maximum_template_length * sizeof(*mapped_positions)
    );
    if (
        projection_cover == NULL
        || insertion_cover == NULL
        || event_cover == NULL
        || observed_n == NULL
        || mapped_positions == NULL
    ) {
        PyErr_NoMemory();
        goto done;
    }

    for (row = 0; row < fragment_count; ++row) {
        const uint32_t template_flat_begin = load_u32_le(&common[3], row);
        const uint32_t template_length =
            load_u32_le(&common[3], row + 1U) - template_flat_begin;
        const uint32_t reference_begin = load_u32_le(&common[1], row);
        const uint32_t reference_end = load_u32_le(&common[2], row);
        const unsigned char haplotype = load_u8(&common[10], row);
        const uint32_t projection_begin = load_u32_le(&details[0], row);
        const uint32_t projection_end = load_u32_le(&details[0], row + 1U);
        const uint32_t event_begin = load_u32_le(&details[1], row);
        const uint32_t event_end = load_u32_le(&details[1], row + 1U);
        const uint32_t site_begin = load_u32_le(&common[5], row);
        const uint32_t site_end = load_u32_le(&common[5], row + 1U);
        const uint32_t original_n_begin = load_u32_le(&details[2], row);
        const uint32_t original_n_end = load_u32_le(&details[2], row + 1U);
        uint32_t previous_template_end = 0;
        uint32_t previous_reference_begin = 0;
        uint32_t previous_reference_end = 0;
        uint32_t previous_variant_index = 0;
        uint32_t previous_n = 0;
        int has_previous_projection = 0;
        int has_previous_event = 0;
        int has_previous_n = 0;

        memset(projection_cover, 0, (size_t)template_length);
        memset(insertion_cover, 0, (size_t)template_length);
        memset(event_cover, 0, (size_t)template_length);
        memset(observed_n, 0, (size_t)template_length);
        memset(mapped_positions, 0xff, (size_t)template_length * sizeof(*mapped_positions));

        for (index = projection_begin; index < projection_end; ++index) {
            const uint32_t template_begin = load_u32_le(&details[3], index);
            const uint32_t template_end = load_u32_le(&details[4], index);
            const uint32_t mapped_begin = load_u32_le(&details[5], index);
            const uint64_t mapped_end_64 =
                (uint64_t)mapped_begin + (template_end - template_begin);
            uint32_t local_offset;
            if (template_begin >= template_end || template_end > template_length) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "projection run has an invalid template interval"
                );
                goto done;
            }
            if (mapped_end_64 > UINT32_MAX) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "projection run reference interval overflows u32"
                );
                goto done;
            }
            if (mapped_begin < reference_begin || mapped_end_64 > reference_end) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "projection run exceeds its reference envelope"
                );
                goto done;
            }
            if (has_previous_projection) {
                if (template_begin < previous_template_end) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "projection runs overlap or are unordered"
                    );
                    goto done;
                }
                if (mapped_begin < previous_reference_end) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "projection references are not increasing"
                    );
                    goto done;
                }
                if (
                    template_begin == previous_template_end
                    && mapped_begin == previous_reference_end
                ) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "projection runs are not maximal"
                    );
                    goto done;
                }
                if (mapped_begin <= previous_reference_begin) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "projection references are not strictly ordered"
                    );
                    goto done;
                }
            }
            for (local_offset = template_begin; local_offset < template_end; ++local_offset) {
                if (projection_cover[local_offset]) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "projection runs cover a template base twice"
                    );
                    goto done;
                }
                projection_cover[local_offset] = 1;
                mapped_positions[local_offset] =
                    mapped_begin + (local_offset - template_begin);
            }
            previous_template_end = template_end;
            previous_reference_begin = mapped_begin;
            previous_reference_end = (uint32_t)mapped_end_64;
            has_previous_projection = 1;
        }

        for (index = event_begin; index < event_end; ++index) {
            const uint32_t variant_index = load_u32_le(&details[6], index);
            const unsigned char kind = load_u8(&details[17], index);
            const unsigned char phased = load_u8(&details[18], index);
            const uint32_t event_reference_begin = load_u32_le(&details[8], index);
            const uint32_t event_reference_end = load_u32_le(&details[9], index);
            const uint32_t event_template_begin = load_u32_le(&details[10], index);
            const uint32_t event_template_end = load_u32_le(&details[11], index);
            const uint32_t ref_begin = load_u32_le(&details[12], index);
            const uint32_t ref_end = load_u32_le(&details[12], index + 1U);
            const uint32_t alt_begin = load_u32_le(&details[13], index);
            const uint32_t alt_end = load_u32_le(&details[13], index + 1U);
            uint32_t reference_span;
            uint32_t template_span;
            uint32_t relative;

            if (variant_index == NO_REFERENCE_POSITION) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event id uses the reserved sentinel"
                );
                goto done;
            }
            if (has_previous_event && variant_index <= previous_variant_index) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event ids must be strictly increasing"
                );
                goto done;
            }
            previous_variant_index = variant_index;
            has_previous_event = 1;
            if (kind < 1 || kind > 3) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event.kind is outside VariantKind"
                );
                goto done;
            }
            if (phased != 255 && phased != haplotype) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event phased haplotype disagrees with fragment"
                );
                goto done;
            }
            if (
                event_reference_begin > event_reference_end
                || event_reference_begin < reference_begin
                || event_reference_end > reference_end
            ) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event reference interval exceeds its fragment"
                );
                goto done;
            }
            if (
                event_template_begin > event_template_end
                || event_template_end > template_length
            ) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event template interval exceeds its fragment"
                );
                goto done;
            }
            reference_span = event_reference_end - event_reference_begin;
            template_span = event_template_end - event_template_begin;
            if (ref_end - ref_begin != reference_span) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event REF bases disagree with its reference span"
                );
                goto done;
            }
            if (alt_end - alt_begin != template_span) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "event ALT bases disagree with its template span"
                );
                goto done;
            }
            if (kind == 1) {
                if (reference_span == 0 || reference_span != template_span) {
                    PyErr_SetString(PyExc_ValueError, "SNV event has an invalid shape");
                    goto done;
                }
            } else if (kind == 2) {
                if (reference_span != 0 || template_span == 0) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "insertion event has an invalid shape"
                    );
                    goto done;
                }
            } else if (reference_span == 0 || template_span != 0) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "deletion event has an invalid shape"
                );
                goto done;
            }

            if (kind == 2) {
                uint32_t previous_mapped = 0;
                uint32_t next_mapped = 0;
                uint32_t cursor;
                int has_previous_mapped = 0;
                int has_next_mapped = 0;
                for (cursor = event_template_begin; cursor > 0; --cursor) {
                    const uint32_t position = mapped_positions[cursor - 1U];
                    if (position != NO_REFERENCE_POSITION) {
                        previous_mapped = position;
                        has_previous_mapped = 1;
                        break;
                    }
                }
                for (cursor = event_template_end; cursor < template_length; ++cursor) {
                    const uint32_t position = mapped_positions[cursor];
                    if (position != NO_REFERENCE_POSITION) {
                        next_mapped = position;
                        has_next_mapped = 1;
                        break;
                    }
                }
                if (
                    (!has_previous_mapped && event_reference_begin != reference_begin)
                    || (has_previous_mapped && previous_mapped >= event_reference_begin)
                    || (!has_next_mapped && event_reference_begin != reference_end)
                    || (has_next_mapped && next_mapped < event_reference_begin)
                ) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "insertion anchor disagrees with projection"
                    );
                    goto done;
                }
            } else if (kind == 3) {
                uint32_t previous_mapped = 0;
                uint32_t next_mapped = 0;
                uint32_t cursor;
                int has_previous_mapped = 0;
                int has_next_mapped = 0;
                int mapped_inside_event = 0;
                for (cursor = event_template_begin; cursor > 0; --cursor) {
                    const uint32_t position = mapped_positions[cursor - 1U];
                    if (position != NO_REFERENCE_POSITION) {
                        previous_mapped = position;
                        has_previous_mapped = 1;
                        break;
                    }
                }
                for (cursor = event_template_begin; cursor < template_length; ++cursor) {
                    const uint32_t position = mapped_positions[cursor];
                    if (position != NO_REFERENCE_POSITION) {
                        next_mapped = position;
                        has_next_mapped = 1;
                        break;
                    }
                }
                for (cursor = 0; cursor < template_length; ++cursor) {
                    const uint32_t position = mapped_positions[cursor];
                    if (
                        position != NO_REFERENCE_POSITION
                        && position >= event_reference_begin
                        && position < event_reference_end
                    ) {
                        mapped_inside_event = 1;
                        break;
                    }
                }
                if (
                    (has_previous_mapped && previous_mapped >= event_reference_begin)
                    || (has_next_mapped && next_mapped < event_reference_end)
                    || mapped_inside_event
                ) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "deletion boundary disagrees with projection"
                    );
                    goto done;
                }
            }

            for (relative = 0; relative < template_span; ++relative) {
                const uint32_t local_offset = event_template_begin + relative;
                if (event_cover[local_offset]) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "event template spans overlap"
                    );
                    goto done;
                }
                event_cover[local_offset] = 1;
                if (
                    load_u8(&common[17], template_flat_begin + local_offset)
                    != load_u8(&details[21], alt_begin + relative)
                ) {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "event ALT bases disagree with the template"
                    );
                    goto done;
                }
                if (kind == 1) {
                    if (
                        mapped_positions[local_offset]
                        != event_reference_begin + relative
                    ) {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "SNV projection disagrees with its event"
                        );
                        goto done;
                    }
                } else if (kind == 2) {
                    if (projection_cover[local_offset]) {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "insertion overlaps a projection run"
                        );
                        goto done;
                    }
                    insertion_cover[local_offset] = 1;
                }
            }
        }

        for (index = 0; index < template_length; ++index) {
            if (projection_cover[index] + insertion_cover[index] != 1) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "projection and insertion runs do not cover the template exactly"
                );
                goto done;
            }
        }
        for (index = site_begin; index < site_end; ++index) {
            const uint32_t local_offset = load_u32_le(&common[8], index);
            if (load_u32_le(&details[14], index) != mapped_positions[local_offset]) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "site reference position disagrees with projection"
                );
                goto done;
            }
        }
        for (index = original_n_begin; index < original_n_end; ++index) {
            const uint32_t local_offset = load_u32_le(&details[15], index);
            if (local_offset >= template_length) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "original-N offset exceeds its template"
                );
                goto done;
            }
            if (has_previous_n && local_offset <= previous_n) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "original-N offsets must be strictly increasing"
                );
                goto done;
            }
            previous_n = local_offset;
            has_previous_n = 1;
            observed_n[local_offset] = 1;
            if (load_u8(&common[17], template_flat_begin + local_offset) != 4) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "PRESERVE_N provenance does not point to N"
                );
                goto done;
            }
        }
        for (index = 0; index < template_length; ++index) {
            const int is_n = load_u8(&common[17], template_flat_begin + index) == 4;
            if ((observed_n[index] != 0) != is_n) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "PRESERVE_N provenance is incomplete"
                );
                goto done;
            }
        }
    }
    result = 1;

done:
    PyMem_Free(projection_cover);
    PyMem_Free(insertion_cover);
    PyMem_Free(event_cover);
    PyMem_Free(observed_n);
    PyMem_Free(mapped_positions);
    return result;
}


static void
store_u32_le(unsigned char *data, uint64_t index, uint32_t value)
{
    const uint64_t offset = index * 4U;
    data[offset] = (unsigned char)(value & 0xffU);
    data[offset + 1U] = (unsigned char)((value >> 8) & 0xffU);
    data[offset + 2U] = (unsigned char)((value >> 16) & 0xffU);
    data[offset + 3U] = (unsigned char)((value >> 24) & 0xffU);
}


static void
store_u64_le(unsigned char *data, uint64_t index, uint64_t value)
{
    const uint64_t offset = index * 8U;
    uint32_t byte_index;
    for (byte_index = 0; byte_index < 8; ++byte_index) {
        data[offset + byte_index] =
            (unsigned char)((value >> (byte_index * 8U)) & 0xffU);
    }
}


static int
checked_byte_length(uint64_t count, uint32_t item_size, Py_ssize_t *result)
{
    const uint64_t length = count * item_size;
    if (length > (uint64_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "packed columns exceed Py_ssize_t");
        return 0;
    }
    *result = (Py_ssize_t)length;
    return 1;
}


static PyObject *
copy_column_bytes(const ColumnBuffer *column)
{
    return PyBytes_FromStringAndSize(
        (const char *)column->view.buf,
        column->view.len
    );
}


static PyObject *
u64_bytes_from_u32_column(const ColumnBuffer *column, uint32_t count)
{
    PyObject *result;
    unsigned char *data;
    Py_ssize_t length;
    uint32_t index;
    if (!checked_byte_length(count, 8, &length)) {
        return NULL;
    }
    result = PyBytes_FromStringAndSize(NULL, length);
    if (result == NULL) {
        return NULL;
    }
    data = (unsigned char *)PyBytes_AS_STRING(result);
    for (index = 0; index < count; ++index) {
        store_u64_le(data, index, load_u32_le(column, index));
    }
    return result;
}


static uint32_t
lower_bound_site_offset(
    const ColumnBuffer *site_offsets,
    uint32_t begin,
    uint32_t end,
    uint32_t target
)
{
    while (begin < end) {
        const uint32_t middle = begin + (end - begin) / 2U;
        if (load_u32_le(site_offsets, middle) < target) {
            begin = middle + 1U;
        } else {
            end = middle;
        }
    }
    return begin;
}


static void
decref_object_array(PyObject **values, Py_ssize_t count)
{
    Py_ssize_t index;
    for (index = 0; index < count; ++index) {
        Py_XDECREF(values[index]);
        values[index] = NULL;
    }
}


static PyObject *
tuple_from_owned_objects(PyObject **values, Py_ssize_t count)
{
    PyObject *result = PyTuple_New(count);
    Py_ssize_t index;
    if (result == NULL) {
        return NULL;
    }
    for (index = 0; index < count; ++index) {
        PyTuple_SET_ITEM(result, index, values[index]);
        values[index] = NULL;
    }
    return result;
}


PyObject *
bsreadsim_cext_pack_protocol_common_columns(PyObject *self, PyObject *args)
{
    PyObject *common_object;
    PyObject *first_ordinal_object;
    PyObject *contig_names_object;
    PyObject *contig_lengths_object;
    PyObject *mates_per_fragment_object;
    PyObject *read_length_r1_object;
    PyObject *read_length_r2_object;
    ColumnBuffer common[COMMON_COLUMN_COUNT];
    uint32_t *contig_lengths = NULL;
    uint32_t contig_count;
    uint32_t first_ordinal;
    uint32_t mates_per_fragment;
    uint32_t read_length_r1;
    uint32_t read_length_r2;
    uint32_t fragment_count;
    uint32_t template_count;
    uint32_t mate_count;
    uint32_t site_count;
    uint32_t maximum_template_length;
    PyObject *ordinals = NULL;
    PyObject *selected_contig_names = NULL;
    PyObject *model_values[12];
    PyObject *fragment_values[12];
    PyObject *model_columns = NULL;
    PyObject *fragment_columns = NULL;
    PyObject *result = NULL;
    Py_ssize_t index;
    Py_ssize_t length;
    uint32_t row;
    uint64_t site_reference_count = 0;
    unsigned char *ordinal_data;
    unsigned char *site_index_data;
    unsigned char *site_reference_offset_data;
    unsigned char *site_reference_read_data;
    unsigned char *site_reference_index_data;
    unsigned char *mate_reference_start_data;
    unsigned char *mate_reference_end_data;
    (void)self;

    memset(common, 0, sizeof(common));
    memset(model_values, 0, sizeof(model_values));
    memset(fragment_values, 0, sizeof(fragment_values));
    if (!PyArg_UnpackTuple(
        args,
        "pack_protocol_common_columns",
        7,
        7,
        &common_object,
        &first_ordinal_object,
        &contig_names_object,
        &contig_lengths_object,
        &mates_per_fragment_object,
        &read_length_r1_object,
        &read_length_r2_object
    )) {
        return NULL;
    }
    if (
        !parse_u32(first_ordinal_object, "batch.first_fragment_ordinal", &first_ordinal)
        || !parse_u32(
            mates_per_fragment_object,
            "header.mates_per_fragment",
            &mates_per_fragment
        )
        || !parse_u32(read_length_r1_object, "header.read_length_r1", &read_length_r1)
        || !parse_u32(read_length_r2_object, "header.read_length_r2", &read_length_r2)
    ) {
        goto done;
    }
    if (!PyTuple_Check(contig_names_object) || !PyTuple_Check(contig_lengths_object)) {
        PyErr_SetString(
            PyExc_TypeError,
            "contig names and lengths must be tuples"
        );
        goto done;
    }
    if (
        PyTuple_GET_SIZE(contig_lengths_object) <= 0
        || PyTuple_GET_SIZE(contig_lengths_object) != PyTuple_GET_SIZE(contig_names_object)
        || (uint64_t)PyTuple_GET_SIZE(contig_lengths_object) > UINT32_MAX
    ) {
        PyErr_SetString(
            PyExc_ValueError,
            "contig names and lengths disagree"
        );
        goto done;
    }
    contig_count = (uint32_t)PyTuple_GET_SIZE(contig_lengths_object);
    contig_lengths = PyMem_Malloc((size_t)contig_count * sizeof(*contig_lengths));
    if (contig_lengths == NULL) {
        PyErr_NoMemory();
        goto done;
    }
    for (index = 0; index < (Py_ssize_t)contig_count; ++index) {
        if (
            !PyUnicode_Check(PyTuple_GET_ITEM(contig_names_object, index))
            || !parse_u32(
                PyTuple_GET_ITEM(contig_lengths_object, index),
                "header.contig.length",
                &contig_lengths[index]
            )
        ) {
            if (!PyErr_Occurred()) {
                PyErr_SetString(PyExc_TypeError, "contig name must be text");
            }
            goto done;
        }
        if (contig_lengths[index] == 0) {
            PyErr_SetString(PyExc_ValueError, "header contig length must be positive");
            goto done;
        }
    }
    if (!acquire_columns(
        common_object,
        COMMON_COLUMN_COUNT,
        common_column_names,
        common,
        "common_columns"
    )) {
        goto done;
    }
    if (!validate_common_columns(
        common,
        contig_lengths,
        contig_count,
        first_ordinal,
        mates_per_fragment,
        read_length_r1,
        read_length_r2,
        0,
        0,
        &fragment_count,
        &template_count,
        &mate_count,
        &site_count,
        &maximum_template_length
    )) {
        goto done;
    }
    (void)template_count;
    (void)maximum_template_length;

    ordinals = PyTuple_New(fragment_count);
    selected_contig_names = PyTuple_New(fragment_count);
    if (ordinals == NULL || selected_contig_names == NULL) {
        goto done;
    }
    if (!checked_byte_length(fragment_count, 8, &length)) {
        goto done;
    }
    model_values[0] = PyBytes_FromStringAndSize(NULL, length);
    if (model_values[0] == NULL) {
        goto done;
    }
    ordinal_data = (unsigned char *)PyBytes_AS_STRING(model_values[0]);
    for (row = 0; row < fragment_count; ++row) {
        const uint32_t ordinal = first_ordinal + row;
        const uint32_t contig_index = load_u32_le(&common[0], row);
        PyObject *ordinal_object = PyLong_FromUnsignedLong(ordinal);
        PyObject *name = PyTuple_GET_ITEM(contig_names_object, contig_index);
        if (ordinal_object == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(ordinals, row, ordinal_object);
        Py_INCREF(name);
        PyTuple_SET_ITEM(selected_contig_names, row, name);
        store_u64_le(ordinal_data, row, ordinal);
    }

    model_values[1] = copy_column_bytes(&common[0]);
    model_values[2] = u64_bytes_from_u32_column(&common[3], fragment_count + 1U);
    model_values[3] = copy_column_bytes(&common[17]);
    model_values[4] = u64_bytes_from_u32_column(&common[5], fragment_count + 1U);
    if (!checked_byte_length(site_count, 4, &length)) {
        goto done;
    }
    model_values[5] = PyBytes_FromStringAndSize(NULL, length);
    model_values[6] = copy_column_bytes(&common[8]);
    if (!checked_byte_length(site_count, 8, &length)) {
        goto done;
    }
    model_values[7] = PyBytes_FromStringAndSize(NULL, length);
    model_values[8] = copy_column_bytes(&common[14]);
    model_values[9] = copy_column_bytes(&common[15]);
    model_values[10] = copy_column_bytes(&common[16]);
    model_values[11] = copy_column_bytes(&common[9]);
    for (index = 1; index < 12; ++index) {
        if (model_values[index] == NULL) {
            goto done;
        }
    }
    site_index_data = (unsigned char *)PyBytes_AS_STRING(model_values[5]);
    memset(PyBytes_AS_STRING(model_values[7]), 0xff, (size_t)length);
    for (row = 0; row < fragment_count; ++row) {
        const uint32_t site_begin = load_u32_le(&common[5], row);
        const uint32_t site_end = load_u32_le(&common[5], row + 1U);
        uint32_t site_index;
        for (site_index = site_begin; site_index < site_end; ++site_index) {
            store_u32_le(site_index_data, site_index, site_index - site_begin);
        }
    }

    fragment_values[0] = copy_column_bytes(&common[10]);
    fragment_values[1] = copy_column_bytes(&common[11]);
    fragment_values[2] = u64_bytes_from_u32_column(&common[4], fragment_count + 1U);
    fragment_values[3] = copy_column_bytes(&common[12]);
    fragment_values[4] = copy_column_bytes(&common[13]);
    fragment_values[5] = copy_column_bytes(&common[6]);
    fragment_values[6] = copy_column_bytes(&common[7]);
    if (!checked_byte_length(mate_count, 8, &length)) {
        goto done;
    }
    fragment_values[7] = PyBytes_FromStringAndSize(NULL, length);
    fragment_values[8] = PyBytes_FromStringAndSize(NULL, length);
    if (!checked_byte_length((uint64_t)mate_count + 1U, 8, &length)) {
        goto done;
    }
    fragment_values[9] = PyBytes_FromStringAndSize(NULL, length);
    for (index = 0; index < 10; ++index) {
        if (fragment_values[index] == NULL) {
            goto done;
        }
    }
    mate_reference_start_data =
        (unsigned char *)PyBytes_AS_STRING(fragment_values[7]);
    mate_reference_end_data =
        (unsigned char *)PyBytes_AS_STRING(fragment_values[8]);
    site_reference_offset_data =
        (unsigned char *)PyBytes_AS_STRING(fragment_values[9]);
    store_u64_le(site_reference_offset_data, 0, 0);
    for (row = 0; row < fragment_count; ++row) {
        const uint32_t reference_begin = load_u32_le(&common[1], row);
        const uint32_t reference_end = load_u32_le(&common[2], row);
        const uint32_t site_begin = load_u32_le(&common[5], row);
        const uint32_t site_end = load_u32_le(&common[5], row + 1U);
        const uint32_t mate_begin = load_u32_le(&common[4], row);
        const uint32_t mate_end = load_u32_le(&common[4], row + 1U);
        uint32_t mate_index;
        for (mate_index = mate_begin; mate_index < mate_end; ++mate_index) {
            const uint32_t template_begin = load_u32_le(&common[6], mate_index);
            const uint32_t template_end = load_u32_le(&common[7], mate_index);
            const uint32_t selected_begin = lower_bound_site_offset(
                &common[8], site_begin, site_end, template_begin
            );
            const uint32_t selected_end = lower_bound_site_offset(
                &common[8], site_begin, site_end, template_end
            );
            site_reference_count += selected_end - selected_begin;
            store_u64_le(
                site_reference_offset_data,
                (uint64_t)mate_index + 1U,
                site_reference_count
            );
            store_u64_le(mate_reference_start_data, mate_index, reference_begin);
            store_u64_le(mate_reference_end_data, mate_index, reference_end);
        }
    }
    if (!checked_byte_length(site_reference_count, 4, &length)) {
        goto done;
    }
    fragment_values[10] = PyBytes_FromStringAndSize(NULL, length);
    fragment_values[11] = PyBytes_FromStringAndSize(NULL, length);
    if (fragment_values[10] == NULL || fragment_values[11] == NULL) {
        goto done;
    }
    site_reference_read_data =
        (unsigned char *)PyBytes_AS_STRING(fragment_values[10]);
    site_reference_index_data =
        (unsigned char *)PyBytes_AS_STRING(fragment_values[11]);
    site_reference_count = 0;
    for (row = 0; row < fragment_count; ++row) {
        const uint32_t site_begin = load_u32_le(&common[5], row);
        const uint32_t site_end = load_u32_le(&common[5], row + 1U);
        const uint32_t mate_begin = load_u32_le(&common[4], row);
        const uint32_t mate_end = load_u32_le(&common[4], row + 1U);
        uint32_t mate_index;
        for (mate_index = mate_begin; mate_index < mate_end; ++mate_index) {
            const uint32_t template_begin = load_u32_le(&common[6], mate_index);
            const uint32_t template_end = load_u32_le(&common[7], mate_index);
            const uint32_t selected_begin = lower_bound_site_offset(
                &common[8], site_begin, site_end, template_begin
            );
            const uint32_t selected_end = lower_bound_site_offset(
                &common[8], site_begin, site_end, template_end
            );
            const int reverse = load_u8(&common[13], mate_index) != 0;
            uint32_t cursor;
            if (reverse) {
                for (cursor = selected_end; cursor > selected_begin; --cursor) {
                    const uint32_t site_index = cursor - 1U;
                    const uint32_t template_offset = load_u32_le(&common[8], site_index);
                    store_u32_le(
                        site_reference_read_data,
                        site_reference_count,
                        template_end - 1U - template_offset
                    );
                    store_u32_le(
                        site_reference_index_data,
                        site_reference_count,
                        site_index - site_begin
                    );
                    site_reference_count += 1U;
                }
            } else {
                for (cursor = selected_begin; cursor < selected_end; ++cursor) {
                    const uint32_t template_offset = load_u32_le(&common[8], cursor);
                    store_u32_le(
                        site_reference_read_data,
                        site_reference_count,
                        template_offset - template_begin
                    );
                    store_u32_le(
                        site_reference_index_data,
                        site_reference_count,
                        cursor - site_begin
                    );
                    site_reference_count += 1U;
                }
            }
        }
    }

    model_columns = tuple_from_owned_objects(model_values, 12);
    fragment_columns = tuple_from_owned_objects(fragment_values, 12);
    if (model_columns == NULL || fragment_columns == NULL) {
        goto done;
    }
    result = PyTuple_Pack(
        4,
        ordinals,
        selected_contig_names,
        model_columns,
        fragment_columns
    );

done:
    release_columns(common, COMMON_COLUMN_COUNT);
    PyMem_Free(contig_lengths);
    decref_object_array(model_values, 12);
    decref_object_array(fragment_values, 12);
    Py_XDECREF(ordinals);
    Py_XDECREF(selected_contig_names);
    Py_XDECREF(model_columns);
    Py_XDECREF(fragment_columns);
    return result;
}


static PyObject *
make_variant(ColumnBuffer *details, uint32_t index)
{
    const uint32_t ref_begin = load_u32_le(&details[12], index);
    const uint32_t ref_end = load_u32_le(&details[12], index + 1U);
    const uint32_t alt_begin = load_u32_le(&details[13], index);
    const uint32_t alt_end = load_u32_le(&details[13], index + 1U);
    const uint32_t id_begin = load_u32_le(&details[7], index);
    const uint32_t id_end = load_u32_le(&details[7], index + 1U);
    const unsigned char source = load_u8(&details[16], index);
    const unsigned char kind = load_u8(&details[17], index);
    PyObject *variant_id = NULL;
    PyObject *ref_bases = NULL;
    PyObject *alt_bases = NULL;
    PyObject *arguments = NULL;
    PyObject *result = NULL;

    ref_bases = PyBytes_FromStringAndSize(
        (const char *)details[20].view.buf + ref_begin,
        (Py_ssize_t)(ref_end - ref_begin)
    );
    alt_bases = PyBytes_FromStringAndSize(
        (const char *)details[21].view.buf + alt_begin,
        (Py_ssize_t)(alt_end - alt_begin)
    );
    variant_id = PyUnicode_DecodeUTF8(
        (const char *)details[19].view.buf + id_begin,
        (Py_ssize_t)(id_end - id_begin),
        "strict"
    );
    arguments = PyTuple_New(9);
    if (
        ref_bases == NULL
        || alt_bases == NULL
        || variant_id == NULL
        || arguments == NULL
        || !tuple_set_unsigned(arguments, 0, load_u32_le(&details[6], index))
        || !tuple_set_borrowed(arguments, 2, variant_sources[source])
        || !tuple_set_borrowed(arguments, 3, variant_kinds[kind])
        || !tuple_set_unsigned(arguments, 4, load_u8(&details[18], index))
        || !tuple_set_unsigned(arguments, 5, load_u32_le(&details[8], index))
        || !tuple_set_unsigned(arguments, 6, load_u32_le(&details[9], index))
    ) {
        goto done;
    }
    PyTuple_SET_ITEM(arguments, 1, variant_id);
    variant_id = NULL;
    PyTuple_SET_ITEM(arguments, 7, ref_bases);
    ref_bases = NULL;
    PyTuple_SET_ITEM(arguments, 8, alt_bases);
    alt_bases = NULL;
    result = call_record(variant_type, arguments);
    arguments = NULL;

done:
    Py_XDECREF(arguments);
    Py_XDECREF(variant_id);
    Py_XDECREF(ref_bases);
    Py_XDECREF(alt_bases);
    return result;
}


static PyObject *
make_methylation_site(
    ColumnBuffer *common,
    ColumnBuffer *details,
    uint32_t global_index,
    uint32_t local_index
)
{
    const uint32_t reference_position = details == NULL
        ? NO_REFERENCE_POSITION
        : load_u32_le(&details[14], global_index);
    const unsigned char context = load_u8(&common[14], global_index);
    const unsigned char source = load_u8(&common[15], global_index);
    const unsigned char allele = load_u8(&common[16], global_index);
    PyObject *arguments = PyTuple_New(7);
    PyObject *probability = PyFloat_FromDouble(
        (double)load_f32_le(&common[9], global_index)
    );
    PyObject *result = NULL;
    if (
        arguments == NULL
        || probability == NULL
        || !tuple_set_unsigned(arguments, 0, local_index)
        || !tuple_set_unsigned(
            arguments,
            1,
            load_u32_le(&common[8], global_index)
        )
        || !tuple_set_signed(
            arguments,
            2,
            reference_position == NO_REFERENCE_POSITION
                ? -1
                : (int64_t)reference_position
        )
        || !tuple_set_borrowed(arguments, 3, methylation_contexts[context])
        || !tuple_set_borrowed(arguments, 4, methylation_sources[source])
        || !tuple_set_borrowed(arguments, 5, methylation_alleles[allele])
    ) {
        goto done;
    }
    PyTuple_SET_ITEM(arguments, 6, probability);
    probability = NULL;
    result = call_record(methylation_site_type, arguments);
    arguments = NULL;

done:
    Py_XDECREF(arguments);
    Py_XDECREF(probability);
    return result;
}


static PyObject *
make_site_reference(uint32_t read_offset, uint32_t site_index)
{
    PyObject *arguments = PyTuple_New(2);
    if (
        arguments == NULL
        || !tuple_set_unsigned(arguments, 0, read_offset)
        || !tuple_set_unsigned(arguments, 1, site_index)
    ) {
        Py_XDECREF(arguments);
        return NULL;
    }
    return call_record(site_reference_type, arguments);
}


static int
find_event_anchor(
    ColumnBuffer *details,
    uint32_t event_begin,
    uint32_t event_end,
    uint32_t variant_index,
    uint32_t *anchor
)
{
    uint32_t index;
    for (index = event_begin; index < event_end; ++index) {
        if (load_u32_le(&details[6], index) == variant_index) {
            *anchor = load_u32_le(&details[8], index);
            return 1;
        }
    }
    PyErr_SetString(PyExc_ValueError, "base event id has no typed event");
    return 0;
}


static PyObject *
make_mate(
    ColumnBuffer *common,
    ColumnBuffer *details,
    uint32_t fragment_reference_start,
    uint32_t fragment_reference_end,
    uint32_t mate_index,
    uint32_t site_begin,
    uint32_t site_end,
    uint32_t event_begin,
    uint32_t event_end,
    const uint32_t *mapped_positions,
    const uint32_t *base_variant_indices
)
{
    const uint32_t template_begin = load_u32_le(&common[6], mate_index);
    const uint32_t template_end = load_u32_le(&common[7], mate_index);
    const int reverse = load_u8(&common[13], mate_index) != 0;
    const uint32_t selected_begin = lower_bound_site_offset(
        &common[8], site_begin, site_end, template_begin
    );
    const uint32_t selected_end = lower_bound_site_offset(
        &common[8], site_begin, site_end, template_end
    );
    const uint32_t site_reference_count = selected_end - selected_begin;
    uint32_t reference_start = details == NULL
        ? fragment_reference_start
        : UINT32_MAX;
    uint32_t reference_end = details == NULL
        ? fragment_reference_end
        : 0;
    int has_mapped_position = details == NULL;
    PyObject *site_references = NULL;
    PyObject *arguments = NULL;
    PyObject *result = NULL;
    uint32_t cursor;
    uint32_t output_index = 0;

    for (cursor = template_begin; cursor < template_end; ++cursor) {
        const uint32_t position = mapped_positions[cursor];
        if (position == NO_REFERENCE_POSITION) {
            continue;
        }
        if (!has_mapped_position || position < reference_start) {
            reference_start = position;
        }
        if (!has_mapped_position || position + 1U > reference_end) {
            reference_end = position + 1U;
        }
        has_mapped_position = 1;
    }
    if (!has_mapped_position) {
        uint32_t anchor = 0;
        int has_anchor = 0;
        for (cursor = template_begin; cursor < template_end; ++cursor) {
            const uint32_t variant_index = base_variant_indices[cursor];
            uint32_t candidate;
            if (variant_index == NO_REFERENCE_POSITION) {
                continue;
            }
            if (!find_event_anchor(
                details,
                event_begin,
                event_end,
                variant_index,
                &candidate
            )) {
                goto done;
            }
            if (has_anchor && candidate != anchor) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "insertion-only v2 mate has ambiguous reference anchors"
                );
                goto done;
            }
            anchor = candidate;
            has_anchor = 1;
        }
        if (!has_anchor) {
            PyErr_SetString(
                PyExc_ValueError,
                "insertion-only v2 mate has no reference anchor"
            );
            goto done;
        }
        reference_start = anchor;
        reference_end = anchor;
    }

    site_references = PyTuple_New(site_reference_count);
    if (site_references == NULL) {
        goto done;
    }
    if (reverse) {
        for (cursor = selected_end; cursor > selected_begin; --cursor) {
            const uint32_t site_index = cursor - 1U;
            const uint32_t template_offset = load_u32_le(&common[8], site_index);
            PyObject *site_reference = make_site_reference(
                template_end - 1U - template_offset,
                site_index - site_begin
            );
            if (site_reference == NULL) {
                goto done;
            }
            PyTuple_SET_ITEM(site_references, output_index, site_reference);
            output_index += 1U;
        }
    } else {
        for (cursor = selected_begin; cursor < selected_end; ++cursor) {
            const uint32_t template_offset = load_u32_le(&common[8], cursor);
            PyObject *site_reference = make_site_reference(
                template_offset - template_begin,
                cursor - site_begin
            );
            if (site_reference == NULL) {
                goto done;
            }
            PyTuple_SET_ITEM(site_references, output_index, site_reference);
            output_index += 1U;
        }
    }
    arguments = PyTuple_New(7);
    if (
        arguments == NULL
        || !tuple_set_unsigned(arguments, 0, load_u8(&common[12], mate_index))
        || !tuple_set_borrowed(arguments, 1, reverse ? Py_True : Py_False)
        || !tuple_set_unsigned(arguments, 2, template_begin)
        || !tuple_set_unsigned(arguments, 3, template_end)
        || !tuple_set_unsigned(arguments, 4, reference_start)
        || !tuple_set_unsigned(arguments, 5, reference_end)
    ) {
        goto done;
    }
    PyTuple_SET_ITEM(arguments, 6, site_references);
    site_references = NULL;
    result = call_record(mate_type, arguments);
    arguments = NULL;

done:
    Py_XDECREF(arguments);
    Py_XDECREF(site_references);
    return result;
}


static PyObject *
make_fragment(
    ColumnBuffer *common,
    ColumnBuffer *details,
    uint32_t first_ordinal,
    uint32_t row,
    uint32_t *mapped_positions,
    uint32_t *base_variant_indices,
    PyObject *no_event_object
)
{
    const uint32_t template_flat_begin = load_u32_le(&common[3], row);
    const uint32_t template_length =
        load_u32_le(&common[3], row + 1U) - template_flat_begin;
    const uint32_t projection_begin =
        details == NULL ? 0U : load_u32_le(&details[0], row);
    const uint32_t projection_end =
        details == NULL ? 0U : load_u32_le(&details[0], row + 1U);
    const uint32_t event_begin =
        details == NULL ? 0U : load_u32_le(&details[1], row);
    const uint32_t event_end =
        details == NULL ? 0U : load_u32_le(&details[1], row + 1U);
    const uint32_t site_begin = load_u32_le(&common[5], row);
    const uint32_t site_end = load_u32_le(&common[5], row + 1U);
    const uint32_t mate_begin = load_u32_le(&common[4], row);
    const uint32_t mate_end = load_u32_le(&common[4], row + 1U);
    PyObject *template_bases = NULL;
    PyObject *reference_positions = NULL;
    PyObject *variant_indices = NULL;
    PyObject *events = NULL;
    PyObject *sites = NULL;
    PyObject *mates = NULL;
    PyObject *arguments = NULL;
    PyObject *result = NULL;
    uint32_t index;

    memset(
        mapped_positions,
        0xff,
        (size_t)template_length * sizeof(*mapped_positions)
    );
    memset(
        base_variant_indices,
        0xff,
        (size_t)template_length * sizeof(*base_variant_indices)
    );
    for (index = projection_begin; index < projection_end; ++index) {
        const uint32_t template_begin = load_u32_le(&details[3], index);
        const uint32_t template_end = load_u32_le(&details[4], index);
        const uint32_t reference_begin = load_u32_le(&details[5], index);
        uint32_t local_offset;
        for (local_offset = template_begin; local_offset < template_end; ++local_offset) {
            mapped_positions[local_offset] =
                reference_begin + (local_offset - template_begin);
        }
    }
    for (index = event_begin; index < event_end; ++index) {
        const uint32_t template_begin = load_u32_le(&details[10], index);
        const uint32_t template_end = load_u32_le(&details[11], index);
        const uint32_t variant_index = load_u32_le(&details[6], index);
        uint32_t local_offset;
        for (local_offset = template_begin; local_offset < template_end; ++local_offset) {
            base_variant_indices[local_offset] = variant_index;
        }
    }

    template_bases = PyBytes_FromStringAndSize(
        (const char *)common[17].view.buf + template_flat_begin,
        template_length
    );
    reference_positions = PyTuple_New(template_length);
    variant_indices = PyTuple_New(template_length);
    events = PyTuple_New(event_end - event_begin);
    sites = PyTuple_New(site_end - site_begin);
    mates = PyTuple_New(mate_end - mate_begin);
    if (
        template_bases == NULL
        || reference_positions == NULL
        || variant_indices == NULL
        || events == NULL
        || sites == NULL
        || mates == NULL
    ) {
        goto done;
    }
    for (index = 0; index < template_length; ++index) {
        const uint32_t position = mapped_positions[index];
        if (!tuple_set_signed(
            reference_positions,
            index,
            position == NO_REFERENCE_POSITION ? -1 : (int64_t)position
        )) {
            goto done;
        }
        if (base_variant_indices[index] == NO_REFERENCE_POSITION) {
            if (!tuple_set_borrowed(variant_indices, index, no_event_object)) {
                goto done;
            }
        } else if (!tuple_set_unsigned(variant_indices, index, base_variant_indices[index])) {
            goto done;
        }
    }
    for (index = event_begin; index < event_end; ++index) {
        PyObject *event = make_variant(details, index);
        if (event == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(events, index - event_begin, event);
    }
    for (index = site_begin; index < site_end; ++index) {
        PyObject *site = make_methylation_site(
            common,
            details,
            index,
            index - site_begin
        );
        if (site == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(sites, index - site_begin, site);
    }
    for (index = mate_begin; index < mate_end; ++index) {
        PyObject *mate = make_mate(
            common,
            details,
            load_u32_le(&common[1], row),
            load_u32_le(&common[2], row),
            index,
            site_begin,
            site_end,
            event_begin,
            event_end,
            mapped_positions,
            base_variant_indices
        );
        if (mate == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(mates, index - mate_begin, mate);
    }

    arguments = PyTuple_New(12);
    if (
        arguments == NULL
        || !tuple_set_unsigned(arguments, 0, first_ordinal + row)
        || !tuple_set_unsigned(arguments, 1, load_u32_le(&common[0], row))
        || !tuple_set_unsigned(arguments, 2, load_u8(&common[10], row))
        || !tuple_set_borrowed(
            arguments,
            3,
            capture_strands[load_u8(&common[11], row)]
        )
        || !tuple_set_unsigned(arguments, 4, load_u32_le(&common[1], row))
        || !tuple_set_unsigned(arguments, 5, load_u32_le(&common[2], row))
    ) {
        goto done;
    }
    PyTuple_SET_ITEM(arguments, 6, template_bases);
    template_bases = NULL;
    PyTuple_SET_ITEM(arguments, 7, reference_positions);
    reference_positions = NULL;
    PyTuple_SET_ITEM(arguments, 8, variant_indices);
    variant_indices = NULL;
    PyTuple_SET_ITEM(arguments, 9, events);
    events = NULL;
    PyTuple_SET_ITEM(arguments, 10, sites);
    sites = NULL;
    PyTuple_SET_ITEM(arguments, 11, mates);
    mates = NULL;
    result = call_record(fragment_type, arguments);
    arguments = NULL;

done:
    Py_XDECREF(arguments);
    Py_XDECREF(mates);
    Py_XDECREF(sites);
    Py_XDECREF(events);
    Py_XDECREF(variant_indices);
    Py_XDECREF(reference_positions);
    Py_XDECREF(template_bases);
    return result;
}


PyObject *
bsreadsim_cext_decode_protocol_fragments(PyObject *self, PyObject *args)
{
    PyObject *common_object;
    PyObject *details_object;
    PyObject *first_ordinal_object;
    PyObject *contig_lengths_object;
    PyObject *mates_per_fragment_object;
    PyObject *read_length_r1_object;
    PyObject *read_length_r2_object;
    ColumnBuffer common[COMMON_COLUMN_COUNT];
    ColumnBuffer details[DETAIL_COLUMN_COUNT];
    uint32_t *contig_lengths = NULL;
    uint32_t *mapped_positions = NULL;
    uint32_t *base_variant_indices = NULL;
    uint32_t contig_count;
    uint32_t first_ordinal;
    uint32_t mates_per_fragment;
    uint32_t read_length_r1;
    uint32_t read_length_r2;
    uint32_t fragment_count;
    uint32_t template_count;
    uint32_t mate_count;
    uint32_t site_count;
    uint32_t maximum_template_length;
    PyObject *no_event_object = NULL;
    PyObject *result = NULL;
    Py_ssize_t index;
    uint32_t row;
    int has_details;
    (void)self;

    memset(common, 0, sizeof(common));
    memset(details, 0, sizeof(details));
    if (!PyArg_UnpackTuple(
        args,
        "decode_protocol_fragments",
        7,
        7,
        &common_object,
        &details_object,
        &first_ordinal_object,
        &contig_lengths_object,
        &mates_per_fragment_object,
        &read_length_r1_object,
        &read_length_r2_object
    )) {
        return NULL;
    }
    has_details = details_object != Py_None;
    if (
        !parse_u32(first_ordinal_object, "batch.first_fragment_ordinal", &first_ordinal)
        || !parse_u32(
            mates_per_fragment_object,
            "header.mates_per_fragment",
            &mates_per_fragment
        )
        || !parse_u32(read_length_r1_object, "header.read_length_r1", &read_length_r1)
        || !parse_u32(read_length_r2_object, "header.read_length_r2", &read_length_r2)
    ) {
        goto done;
    }
    if (!PyTuple_Check(contig_lengths_object)) {
        PyErr_SetString(PyExc_TypeError, "contig_lengths must be a tuple");
        goto done;
    }
    if (
        PyTuple_GET_SIZE(contig_lengths_object) <= 0
        || (uint64_t)PyTuple_GET_SIZE(contig_lengths_object) > UINT32_MAX
    ) {
        PyErr_SetString(PyExc_ValueError, "contig_lengths must not be empty");
        goto done;
    }
    contig_count = (uint32_t)PyTuple_GET_SIZE(contig_lengths_object);
    contig_lengths = PyMem_Malloc((size_t)contig_count * sizeof(*contig_lengths));
    if (contig_lengths == NULL) {
        PyErr_NoMemory();
        goto done;
    }
    for (index = 0; index < (Py_ssize_t)contig_count; ++index) {
        if (!parse_u32(
            PyTuple_GET_ITEM(contig_lengths_object, index),
            "header.contig.length",
            &contig_lengths[index]
        )) {
            goto done;
        }
        if (contig_lengths[index] == 0) {
            PyErr_SetString(PyExc_ValueError, "header contig length must be positive");
            goto done;
        }
    }
    if (!acquire_columns(
        common_object,
        COMMON_COLUMN_COUNT,
        common_column_names,
        common,
        "common_columns"
    )) {
        goto done;
    }
    if (has_details && !acquire_columns(
        details_object,
        DETAIL_COLUMN_COUNT,
        detail_column_names,
        details,
        "has_details"
    )) {
        goto done;
    }
    if (!validate_common_columns(
        common,
        contig_lengths,
        contig_count,
        first_ordinal,
        mates_per_fragment,
        read_length_r1,
        read_length_r2,
        0,
        0,
        &fragment_count,
        &template_count,
        &mate_count,
        &site_count,
        &maximum_template_length
    )) {
        goto done;
    }
    if (has_details && !validate_details(
        common,
        details,
        fragment_count,
        site_count,
        maximum_template_length
    )) {
        goto done;
    }
    (void)template_count;
    (void)mate_count;
    if (!initialize_protocol_types()) {
        goto done;
    }
    mapped_positions = PyMem_Malloc(
        (size_t)maximum_template_length * sizeof(*mapped_positions)
    );
    base_variant_indices = PyMem_Malloc(
        (size_t)maximum_template_length * sizeof(*base_variant_indices)
    );
    no_event_object = PyLong_FromUnsignedLong(NO_REFERENCE_POSITION);
    result = PyTuple_New(fragment_count);
    if (
        mapped_positions == NULL
        || base_variant_indices == NULL
        || no_event_object == NULL
        || result == NULL
    ) {
        if (!PyErr_Occurred()) {
            PyErr_NoMemory();
        }
        Py_CLEAR(result);
        goto done;
    }
    for (row = 0; row < fragment_count; ++row) {
        PyObject *fragment = make_fragment(
            common,
            has_details ? details : NULL,
            first_ordinal,
            row,
            mapped_positions,
            base_variant_indices,
            no_event_object
        );
        if (fragment == NULL) {
            Py_CLEAR(result);
            goto done;
        }
        PyTuple_SET_ITEM(result, row, fragment);
    }

done:
    Py_XDECREF(no_event_object);
    PyMem_Free(base_variant_indices);
    PyMem_Free(mapped_positions);
    PyMem_Free(contig_lengths);
    release_columns(details, DETAIL_COLUMN_COUNT);
    release_columns(common, COMMON_COLUMN_COUNT);
    return result;
}


PyObject *
bsreadsim_cext_validate_protocol_batch_columns(PyObject *self, PyObject *args)
{
    PyObject *common_object;
    PyObject *details_object;
    PyObject *first_ordinal_object;
    PyObject *contig_lengths_object;
    PyObject *mates_per_fragment_object;
    PyObject *read_length_r1_object;
    PyObject *read_length_r2_object;
    PyObject *expected_ordinal_object;
    PyObject *has_details_object;
    ColumnBuffer common[COMMON_COLUMN_COUNT];
    ColumnBuffer details[DETAIL_COLUMN_COUNT];
    uint32_t *contig_lengths = NULL;
    uint32_t contig_count;
    uint32_t first_ordinal;
    uint32_t mates_per_fragment;
    uint32_t read_length_r1;
    uint32_t read_length_r2;
    uint32_t expected_ordinal = 0;
    uint32_t has_details_mode;
    uint32_t fragment_count;
    uint32_t template_count;
    uint32_t mate_count;
    uint32_t site_count;
    uint32_t maximum_template_length;
    Py_ssize_t index;
    int has_expected_ordinal;
    int result = 0;
    (void)self;

    memset(common, 0, sizeof(common));
    memset(details, 0, sizeof(details));
    if (!PyArg_UnpackTuple(
        args,
        "validate_protocol_batch_columns",
        9,
        9,
        &common_object,
        &details_object,
        &first_ordinal_object,
        &contig_lengths_object,
        &mates_per_fragment_object,
        &read_length_r1_object,
        &read_length_r2_object,
        &expected_ordinal_object,
        &has_details_object
    )) {
        return NULL;
    }
    if (
        !parse_u32(first_ordinal_object, "batch.first_fragment_ordinal", &first_ordinal)
        || !parse_u32(
            mates_per_fragment_object,
            "header.mates_per_fragment",
            &mates_per_fragment
        )
        || !parse_u32(read_length_r1_object, "header.read_length_r1", &read_length_r1)
        || !parse_u32(read_length_r2_object, "header.read_length_r2", &read_length_r2)
        || !parse_u32(has_details_object, "header.has_details", &has_details_mode)
    ) {
        goto done;
    }
    if (mates_per_fragment != 1 && mates_per_fragment != 2) {
        PyErr_SetString(
            PyExc_ValueError,
            "header.mates_per_fragment must be one or two"
        );
        goto done;
    }
    if (
        read_length_r1 == 0
        || (mates_per_fragment == 1 && read_length_r2 != 0)
        || (mates_per_fragment == 2 && read_length_r2 == 0)
    ) {
        PyErr_SetString(
            PyExc_ValueError,
            "header read lengths disagree with its mate count"
        );
        goto done;
    }
    if (has_details_mode > 1) {
        PyErr_SetString(PyExc_ValueError, "header.has_details is invalid");
        goto done;
    }
    if (
        (has_details_mode == 0 && details_object != Py_None)
        || (has_details_mode == 1 && details_object == Py_None)
    ) {
        PyErr_SetString(
            PyExc_ValueError,
            "batch details columns disagree with header"
        );
        goto done;
    }

    has_expected_ordinal = expected_ordinal_object != Py_None;
    if (
        has_expected_ordinal
        && !parse_u32(
            expected_ordinal_object,
            "expected_first_ordinal",
            &expected_ordinal
        )
    ) {
        goto done;
    }
    if (!PyTuple_Check(contig_lengths_object)) {
        PyErr_SetString(PyExc_TypeError, "contig_lengths must be a tuple");
        goto done;
    }
    if (
        PyTuple_GET_SIZE(contig_lengths_object) <= 0
        || (uint64_t)PyTuple_GET_SIZE(contig_lengths_object) > UINT32_MAX
    ) {
        PyErr_SetString(PyExc_ValueError, "contig_lengths must not be empty");
        goto done;
    }
    contig_count = (uint32_t)PyTuple_GET_SIZE(contig_lengths_object);
    contig_lengths = PyMem_Malloc((size_t)contig_count * sizeof(*contig_lengths));
    if (contig_lengths == NULL) {
        PyErr_NoMemory();
        goto done;
    }
    for (index = 0; index < (Py_ssize_t)contig_count; ++index) {
        if (!parse_u32(
            PyTuple_GET_ITEM(contig_lengths_object, index),
            "header.contig.length",
            &contig_lengths[index]
        )) {
            goto done;
        }
        if (contig_lengths[index] == 0) {
            PyErr_SetString(PyExc_ValueError, "header contig length must be positive");
            goto done;
        }
    }
    if (!acquire_columns(
        common_object,
        COMMON_COLUMN_COUNT,
        common_column_names,
        common,
        "common_columns"
    )) {
        goto done;
    }
    if (
        has_details_mode == 1
        && !acquire_columns(
            details_object,
            DETAIL_COLUMN_COUNT,
            detail_column_names,
            details,
            "has_details"
        )
    ) {
        goto done;
    }
    if (!validate_common_columns(
        common,
        contig_lengths,
        contig_count,
        first_ordinal,
        mates_per_fragment,
        read_length_r1,
        read_length_r2,
        has_expected_ordinal,
        expected_ordinal,
        &fragment_count,
        &template_count,
        &mate_count,
        &site_count,
        &maximum_template_length
    )) {
        goto done;
    }
    if (
        has_details_mode == 1
        && !validate_details(
            common,
            details,
            fragment_count,
            site_count,
            maximum_template_length
        )
    ) {
        goto done;
    }
    result = 1;

done:
    PyMem_Free(contig_lengths);
    release_columns(details, DETAIL_COLUMN_COUNT);
    release_columns(common, COMMON_COLUMN_COUNT);
    if (!result) {
        return NULL;
    }
    Py_RETURN_NONE;
}
