#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "_native_protocol_types.h"


#define UINT32_MASK UINT32_C(0xffffffff)
#define PHILOX_M0 UINT32_C(0xd2511f53)
#define PHILOX_M1 UINT32_C(0xcd9e8d57)
#define PHILOX_W0 UINT32_C(0x9e3779b9)
#define PHILOX_W1 UINT32_C(0xbb67ae85)


static uint32_t crc32c_tables[8][256];


PyObject *fragment_type;
PyObject *variant_event_type;
PyObject *methylation_site_type;
PyObject *site_reference_type;
PyObject *mate_type;
PyObject *capture_strands[3];
PyObject *variant_kinds[4];
PyObject *methylation_contexts[16];
PyObject *methylation_sources[5];
PyObject *methylation_alleles[3];
static int protocol_types_initialized;
static PyObject *column_site_index_name;
static PyObject *column_probability_name;


static const unsigned char alternative_bases[4][3] = {
    {1, 2, 3},
    {0, 2, 3},
    {0, 1, 3},
    {0, 1, 2}
};


PyObject *bsreadsim_native_truth_json_bytes(PyObject *self, PyObject *args);
PyObject *bsreadsim_native_validate_protocol_batch_columns(
    PyObject *self,
    PyObject *args
);
PyObject *bsreadsim_native_pack_protocol_common_columns(
    PyObject *self,
    PyObject *args
);
PyObject *bsreadsim_native_decode_protocol_fragments(
    PyObject *self,
    PyObject *args
);


static void
initialize_crc32c_table(void)
{
    const uint32_t polynomial = UINT32_C(0x82f63b78);
    uint32_t value;

    for (value = 0; value < 256; ++value) {
        uint32_t crc = value;
        int bit;
        for (bit = 0; bit < 8; ++bit) {
            crc = (crc >> 1) ^ ((crc & 1) ? polynomial : 0);
        }
        crc32c_tables[0][value] = crc;
    }
    for (value = 0; value < 256; ++value) {
        int table;
        uint32_t crc = crc32c_tables[0][value];
        for (table = 1; table < 8; ++table) {
            crc = crc32c_tables[0][crc & UINT32_C(0xff)] ^ (crc >> 8);
            crc32c_tables[table][value] = crc;
        }
    }
}


static void
philox4x32_10(
    uint64_t key,
    uint64_t entity_ordinal,
    uint64_t local_index,
    uint32_t output[4]
)
{
    uint32_t counter_0 = (uint32_t)(entity_ordinal & UINT32_MASK);
    uint32_t counter_1 = (uint32_t)(entity_ordinal >> 32);
    uint32_t counter_2 = (uint32_t)(local_index & UINT32_MASK);
    uint32_t counter_3 = (uint32_t)(local_index >> 32);
    uint32_t key_0 = (uint32_t)(key & UINT32_MASK);
    uint32_t key_1 = (uint32_t)(key >> 32);
    int round;

    for (round = 0; round < 10; ++round) {
        const uint64_t product_0 = (uint64_t)PHILOX_M0 * counter_0;
        const uint64_t product_1 = (uint64_t)PHILOX_M1 * counter_2;
        const uint32_t next_0 =
            (uint32_t)(((product_1 >> 32) ^ counter_1 ^ key_0) & UINT32_MASK);
        const uint32_t next_1 = (uint32_t)(product_1 & UINT32_MASK);
        const uint32_t next_2 =
            (uint32_t)(((product_0 >> 32) ^ counter_3 ^ key_1) & UINT32_MASK);
        const uint32_t next_3 = (uint32_t)(product_0 & UINT32_MASK);

        counter_0 = next_0;
        counter_1 = next_1;
        counter_2 = next_2;
        counter_3 = next_3;
        key_0 = (uint32_t)((key_0 + PHILOX_W0) & UINT32_MASK);
        key_1 = (uint32_t)((key_1 + PHILOX_W1) & UINT32_MASK);
    }

    output[0] = counter_0;
    output[1] = counter_1;
    output[2] = counter_2;
    output[3] = counter_3;
}


static uint64_t
philox_u64(
    uint64_t key,
    uint64_t entity_ordinal,
    uint64_t local_index,
    int pair
)
{
    uint32_t block[4];
    const int offset = pair * 2;
    philox4x32_10(key, entity_ordinal, local_index, block);
    return (uint64_t)block[offset] | ((uint64_t)block[offset + 1] << 32);
}


static int
parse_pair(int pair)
{
    if (pair != 0 && pair != 1) {
        PyErr_SetString(PyExc_ValueError, "pair must be 0 or 1");
        return 0;
    }
    return 1;
}


static PyObject *
native_u64(PyObject *self, PyObject *args)
{
    unsigned long long key;
    unsigned long long entity_ordinal;
    unsigned long long local_index;
    int pair = 0;
    uint64_t value;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "KKK|i:u64",
            &key,
            &entity_ordinal,
            &local_index,
            &pair
        )) {
        return NULL;
    }
    if (!parse_pair(pair)) {
        return NULL;
    }
    value = philox_u64(
        (uint64_t)key,
        (uint64_t)entity_ordinal,
        (uint64_t)local_index,
        pair
    );
    return PyLong_FromUnsignedLongLong((unsigned long long)value);
}


static PyObject *
native_bernoulli(PyObject *self, PyObject *args)
{
    unsigned long long key;
    unsigned long long entity_ordinal;
    unsigned long long local_index;
    double probability;
    int pair = 0;
    uint64_t random_value;
    double uniform;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "KKKd|i:bernoulli",
            &key,
            &entity_ordinal,
            &local_index,
            &probability,
            &pair
        )) {
        return NULL;
    }
    if (!parse_pair(pair)) {
        return NULL;
    }
    if (probability == 0.0) {
        Py_RETURN_FALSE;
    }
    if (probability == 1.0) {
        Py_RETURN_TRUE;
    }

    random_value = philox_u64(
        (uint64_t)key,
        (uint64_t)entity_ordinal,
        (uint64_t)local_index,
        pair
    );
    uniform = ldexp((double)(random_value >> 11), -53);
    if (uniform < probability) {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}


static PyObject *
native_philox_pairs(PyObject *self, PyObject *args)
{
    unsigned long long key;
    Py_buffer entity_ordinals = {0};
    Py_buffer local_indices = {0};
    PyObject *pair_0 = NULL;
    PyObject *pair_1 = NULL;
    PyObject *result = NULL;
    unsigned char *pair_0_data;
    unsigned char *pair_1_data;
    Py_ssize_t offset;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "Ky*y*:philox_pairs",
            &key,
            &entity_ordinals,
            &local_indices
        )) {
        return NULL;
    }
    if (entity_ordinals.len != local_indices.len) {
        PyErr_SetString(
            PyExc_ValueError,
            "Philox counter buffers must have the same length"
        );
        goto done;
    }
    if (entity_ordinals.len % (Py_ssize_t)sizeof(uint64_t) != 0) {
        PyErr_SetString(
            PyExc_ValueError,
            "Philox counter buffers must contain packed uint64 values"
        );
        goto done;
    }

    pair_0 = PyBytes_FromStringAndSize(NULL, entity_ordinals.len);
    pair_1 = PyBytes_FromStringAndSize(NULL, entity_ordinals.len);
    if (pair_0 == NULL || pair_1 == NULL) {
        goto done;
    }
    pair_0_data = (unsigned char *)PyBytes_AS_STRING(pair_0);
    pair_1_data = (unsigned char *)PyBytes_AS_STRING(pair_1);

    Py_BEGIN_ALLOW_THREADS
    for (offset = 0; offset < entity_ordinals.len; offset += sizeof(uint64_t)) {
        uint64_t entity_ordinal;
        uint64_t local_index;
        uint64_t first;
        uint64_t second;
        uint32_t block[4];
        memcpy(
            &entity_ordinal,
            (const unsigned char *)entity_ordinals.buf + offset,
            sizeof(entity_ordinal)
        );
        memcpy(
            &local_index,
            (const unsigned char *)local_indices.buf + offset,
            sizeof(local_index)
        );
        philox4x32_10(
            (uint64_t)key,
            entity_ordinal,
            local_index,
            block
        );
        first = (uint64_t)block[0] | ((uint64_t)block[1] << 32);
        second = (uint64_t)block[2] | ((uint64_t)block[3] << 32);
        memcpy(pair_0_data + offset, &first, sizeof(first));
        memcpy(pair_1_data + offset, &second, sizeof(second));
    }
    Py_END_ALLOW_THREADS

    result = PyTuple_Pack(2, pair_0, pair_1);

done:
    Py_XDECREF(pair_1);
    Py_XDECREF(pair_0);
    if (local_indices.obj != NULL) {
        PyBuffer_Release(&local_indices);
    }
    if (entity_ordinals.obj != NULL) {
        PyBuffer_Release(&entity_ordinals);
    }
    return result;
}


static unsigned int
multiply_high_by_three(uint64_t value)
{
    const uint64_t doubled = value + value;
    const unsigned int first_carry = doubled < value ? 1U : 0U;
    const uint64_t tripled = doubled + value;
    const unsigned int second_carry = tripled < doubled ? 1U : 0U;
    return first_carry + second_carry;
}


static PyObject *
native_apply_uniform_errors(PyObject *self, PyObject *args)
{
    Py_buffer input;
    unsigned long long key;
    unsigned long long entity_ordinal;
    unsigned long long mate_index;
    double probability;
    PyObject *bases = NULL;
    PyObject *flags = NULL;
    PyObject *result = NULL;
    unsigned char *base_data;
    unsigned char *flag_data;
    Py_ssize_t offset;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "y*KKKd:apply_uniform_errors",
            &input,
            &key,
            &entity_ordinal,
            &mate_index,
            &probability
        )) {
        return NULL;
    }
    if (mate_index > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError, "mate_index exceeds uint32");
        goto done;
    }
    if (input.len > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError, "read length exceeds uint32");
        goto done;
    }
    if (!isfinite(probability) || probability < 0.0 || probability > 1.0) {
        PyErr_SetString(
            PyExc_ValueError,
            "error probability must be finite and in [0, 1]"
        );
        goto done;
    }

    bases = PyBytes_FromStringAndSize(NULL, input.len);
    flags = PyBytes_FromStringAndSize(NULL, input.len);
    if (bases == NULL || flags == NULL) {
        goto done;
    }
    base_data = (unsigned char *)PyBytes_AS_STRING(bases);
    flag_data = (unsigned char *)PyBytes_AS_STRING(flags);
    memcpy(base_data, input.buf, (size_t)input.len);
    memset(flag_data, 0, (size_t)input.len);

    for (offset = 0; offset < input.len; ++offset) {
        const unsigned char base = base_data[offset];
        uint32_t block[4];
        uint64_t decision;
        uint64_t alternative;
        uint64_t local_index;
        double uniform;
        int has_error;

        if (base > 4) {
            PyErr_SetString(PyExc_ValueError, "base code must be in [0, 4]");
            goto done;
        }
        if (base == 4 || probability == 0.0) {
            continue;
        }
        local_index = ((uint64_t)mate_index << 32) | (uint32_t)offset;
        philox4x32_10(
            (uint64_t)key,
            (uint64_t)entity_ordinal,
            local_index,
            block
        );
        decision = (uint64_t)block[0] | ((uint64_t)block[1] << 32);
        uniform = ldexp((double)(decision >> 11), -53);
        has_error = probability == 1.0 || uniform < probability;
        if (!has_error) {
            continue;
        }
        alternative = (uint64_t)block[2] | ((uint64_t)block[3] << 32);
        base_data[offset] = alternative_bases[base][
            multiply_high_by_three(alternative)
        ];
        flag_data[offset] = 1;
    }
    result = PyTuple_Pack(2, bases, flags);

done:
    Py_XDECREF(flags);
    Py_XDECREF(bases);
    PyBuffer_Release(&input);
    return result;
}


static uint32_t
crc32c_bytes(const unsigned char *data, Py_ssize_t length)
{
    uint32_t crc = UINT32_C(0xffffffff);
    while (length >= 8) {
        const uint32_t first = crc
            ^ (uint32_t)data[0]
            ^ ((uint32_t)data[1] << 8)
            ^ ((uint32_t)data[2] << 16)
            ^ ((uint32_t)data[3] << 24);
        const uint32_t second =
            (uint32_t)data[4]
            | ((uint32_t)data[5] << 8)
            | ((uint32_t)data[6] << 16)
            | ((uint32_t)data[7] << 24);
        crc =
            crc32c_tables[7][first & UINT32_C(0xff)]
            ^ crc32c_tables[6][(first >> 8) & UINT32_C(0xff)]
            ^ crc32c_tables[5][(first >> 16) & UINT32_C(0xff)]
            ^ crc32c_tables[4][first >> 24]
            ^ crc32c_tables[3][second & UINT32_C(0xff)]
            ^ crc32c_tables[2][(second >> 8) & UINT32_C(0xff)]
            ^ crc32c_tables[1][(second >> 16) & UINT32_C(0xff)]
            ^ crc32c_tables[0][second >> 24];
        data += 8;
        length -= 8;
    }
    while (length != 0) {
        crc = crc32c_tables[0][(crc ^ *data) & UINT32_C(0xff)] ^ (crc >> 8);
        ++data;
        --length;
    }
    return crc ^ UINT32_C(0xffffffff);
}


static PyObject *
native_crc32c(PyObject *self, PyObject *args)
{
    Py_buffer buffer;
    uint32_t value;

    (void)self;
    if (!PyArg_ParseTuple(args, "y*:crc32c", &buffer)) {
        return NULL;
    }
    if (buffer.len >= 4096) {
        Py_BEGIN_ALLOW_THREADS
        value = crc32c_bytes((const unsigned char *)buffer.buf, buffer.len);
        Py_END_ALLOW_THREADS
    } else {
        value = crc32c_bytes((const unsigned char *)buffer.buf, buffer.len);
    }
    PyBuffer_Release(&buffer);
    return PyLong_FromUnsignedLong((unsigned long)value);
}


static PyObject *
get_required_attribute(PyObject *owner, const char *name)
{
    PyObject *value = PyObject_GetAttrString(owner, name);
    if (value == NULL) {
        PyErr_Format(
            PyExc_RuntimeError,
            "native adapter could not resolve bsreadsim.model.%s",
            name
        );
    }
    return value;
}


static PyObject *
make_enum_value(PyObject *enum_type, long value)
{
    PyObject *argument = PyLong_FromLong(value);
    PyObject *result;
    if (argument == NULL) {
        return NULL;
    }
    result = PyObject_CallFunctionObjArgs(enum_type, argument, NULL);
    Py_DECREF(argument);
    return result;
}


int
initialize_protocol_types(void)
{
    PyObject *module = NULL;
    PyObject *capture_type = NULL;
    PyObject *variant_kind_type = NULL;
    PyObject *context_type = NULL;
    PyObject *source_type = NULL;
    PyObject *allele_type = NULL;
    const int context_values[] = {1, 3, 7, 9, 11, 15};
    int index;

    if (protocol_types_initialized) {
        return 1;
    }
    module = PyImport_ImportModule("bsreadsim.model");
    if (module == NULL) {
        return 0;
    }
    fragment_type = get_required_attribute(module, "Fragment");
    variant_event_type = get_required_attribute(module, "VariantEvent");
    methylation_site_type = get_required_attribute(module, "MethylationSite");
    site_reference_type = get_required_attribute(module, "SiteReference");
    mate_type = get_required_attribute(module, "Mate");
    capture_type = get_required_attribute(module, "CaptureStrand");
    variant_kind_type = get_required_attribute(module, "VariantKind");
    context_type = get_required_attribute(module, "MethylationContext");
    source_type = get_required_attribute(module, "MethylationSource");
    allele_type = get_required_attribute(module, "MethylationAllele");
    if (
        fragment_type == NULL
        || variant_event_type == NULL
        || methylation_site_type == NULL
        || site_reference_type == NULL
        || mate_type == NULL
        || capture_type == NULL
        || variant_kind_type == NULL
        || context_type == NULL
        || source_type == NULL
        || allele_type == NULL
    ) {
        goto fail;
    }
    for (index = 0; index < 3; ++index) {
        capture_strands[index] = make_enum_value(capture_type, index);
        if (capture_strands[index] == NULL) {
            goto fail;
        }
    }
    for (index = 1; index <= 3; ++index) {
        variant_kinds[index] = make_enum_value(variant_kind_type, index);
        if (variant_kinds[index] == NULL) {
            goto fail;
        }
    }
    for (index = 0; index < 6; ++index) {
        const int value = context_values[index];
        methylation_contexts[value] = make_enum_value(context_type, value);
        if (methylation_contexts[value] == NULL) {
            goto fail;
        }
    }
    for (index = 1; index <= 4; ++index) {
        methylation_sources[index] = make_enum_value(source_type, index);
        if (methylation_sources[index] == NULL) {
            goto fail;
        }
    }
    for (index = 0; index <= 2; ++index) {
        methylation_alleles[index] = make_enum_value(allele_type, index);
        if (methylation_alleles[index] == NULL) {
            goto fail;
        }
    }
    Py_DECREF(allele_type);
    Py_DECREF(source_type);
    Py_DECREF(context_type);
    Py_DECREF(variant_kind_type);
    Py_DECREF(capture_type);
    Py_DECREF(module);
    protocol_types_initialized = 1;
    return 1;

fail:
    Py_XDECREF(allele_type);
    Py_XDECREF(source_type);
    Py_XDECREF(context_type);
    Py_XDECREF(variant_kind_type);
    Py_XDECREF(capture_type);
    Py_DECREF(module);
    return 0;
}


int
tuple_set_unsigned(PyObject *tuple, Py_ssize_t index, uint64_t value)
{
    PyObject *item = PyLong_FromUnsignedLongLong((unsigned long long)value);
    if (item == NULL) {
        return 0;
    }
    PyTuple_SET_ITEM(tuple, index, item);
    return 1;
}


int
tuple_set_signed(PyObject *tuple, Py_ssize_t index, int64_t value)
{
    PyObject *item = PyLong_FromLongLong((long long)value);
    if (item == NULL) {
        return 0;
    }
    PyTuple_SET_ITEM(tuple, index, item);
    return 1;
}


int
tuple_set_borrowed(PyObject *tuple, Py_ssize_t index, PyObject *value)
{
    Py_INCREF(value);
    PyTuple_SET_ITEM(tuple, index, value);
    return 1;
}


PyObject *
call_record(PyObject *record_type, PyObject *arguments)
{
    PyObject *result = PyObject_CallObject(record_type, arguments);
    Py_DECREF(arguments);
    return result;
}


static int
checked_add_size(Py_ssize_t *total, Py_ssize_t amount, const char *name)
{
    if (amount < 0 || *total > PY_SSIZE_T_MAX - amount) {
        PyErr_Format(PyExc_OverflowError, "%s exceeds the native buffer limit", name);
        return 0;
    }
    *total += amount;
    return 1;
}


static int
checked_buffer_size(Py_ssize_t count, Py_ssize_t width, Py_ssize_t *result)
{
    if (count < 0 || width < 0 || (width != 0 && count > PY_SSIZE_T_MAX / width)) {
        PyErr_SetString(PyExc_OverflowError, "batch exceeds the native buffer limit");
        return 0;
    }
    *result = count * width;
    return 1;
}


static uint64_t
load_u64_le_bytes(const unsigned char *data)
{
    uint64_t value = 0;
    int byte;
    for (byte = 7; byte >= 0; --byte) {
        value = (value << 8) | data[byte];
    }
    return value;
}


static Py_ssize_t
u64_decimal_length(uint64_t value)
{
    Py_ssize_t length = 1;
    while (value >= UINT64_C(10)) {
        value /= UINT64_C(10);
        length += 1;
    }
    return length;
}


static unsigned char *
write_u64_decimal(unsigned char *target, uint64_t value)
{
    unsigned char reversed[20];
    Py_ssize_t length = 0;
    Py_ssize_t index;
    do {
        reversed[length++] = (unsigned char)('0' + value % UINT64_C(10));
        value /= UINT64_C(10);
    } while (value != 0);
    for (index = length; index > 0; --index) {
        *target++ = reversed[index - 1];
    }
    return target;
}


static int
valid_fastq_contig_name(
    PyObject *value,
    const char **utf8,
    Py_ssize_t *utf8_length
)
{
    Py_ssize_t length;
    Py_ssize_t index;
    if (!PyUnicode_Check(value)) {
        return 0;
    }
    length = PyUnicode_GetLength(value);
    if (length <= 0) {
        return length < 0 ? -1 : 0;
    }
    for (index = 0; index < length; ++index) {
        const Py_UCS4 character = PyUnicode_ReadChar(value, index);
        if (character == (Py_UCS4)-1 && PyErr_Occurred()) {
            return -1;
        }
        if (!Py_UNICODE_ISPRINTABLE(character)
            || Py_UNICODE_ISSPACE(character)) {
            return 0;
        }
    }
    *utf8 = PyUnicode_AsUTF8AndSize(value, utf8_length);
    return *utf8 == NULL ? -1 : 1;
}


static unsigned char *
write_fastq_record(
    unsigned char *target,
    const char *contig_name,
    Py_ssize_t contig_name_length,
    uint64_t reference_start,
    uint64_t reference_end,
    uint64_t ordinal,
    unsigned int mate_number,
    const unsigned char *sequence,
    Py_ssize_t read_length,
    unsigned char quality
)
{
    const uint64_t left = reference_start + UINT64_C(1);
    const uint64_t right = reference_end > reference_start
        ? reference_end
        : left;
    *target++ = '@';
    memcpy(target, contig_name, (size_t)contig_name_length);
    target += contig_name_length;
    *target++ = ':';
    target = write_u64_decimal(target, left);
    *target++ = '-';
    target = write_u64_decimal(target, right);
    *target++ = ':';
    target = write_u64_decimal(target, ordinal);
    *target++ = '/';
    *target++ = (unsigned char)('0' + mate_number);
    *target++ = '\n';
    memcpy(target, sequence, (size_t)read_length);
    target += read_length;
    *target++ = '\n';
    *target++ = '+';
    *target++ = '\n';
    memset(target, quality, (size_t)read_length);
    target += read_length;
    *target++ = '\n';
    return target;
}


static PyObject *
native_format_fastq_batch(PyObject *self, PyObject *args)
{
    PyObject *contig_names;
    Py_buffer ordinals = {0};
    Py_buffer mate_offsets = {0};
    Py_buffer mate_indices = {0};
    Py_buffer reference_starts = {0};
    Py_buffer reference_ends = {0};
    Py_buffer sequences = {0};
    Py_ssize_t read_length;
    int quality;
    int paired_end;
    Py_ssize_t fragment_count;
    Py_ssize_t mates_per_fragment;
    Py_ssize_t mate_count;
    Py_ssize_t expected_size;
    Py_ssize_t record_base_size = 0;
    Py_ssize_t read1_size = 0;
    Py_ssize_t read2_size = 0;
    Py_ssize_t fragment_index;
    PyObject *read1 = NULL;
    PyObject *read2 = NULL;
    PyObject *record_lengths = NULL;
    PyObject *result = NULL;
    unsigned char *read1_cursor;
    unsigned char *read2_cursor = NULL;
    const unsigned char *ordinal_data;
    const unsigned char *offset_data;
    const unsigned char *mate_index_data;
    const unsigned char *reference_start_data;
    const unsigned char *reference_end_data;
    const unsigned char *sequence_data;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "Oy*y*y*y*y*y*nii:format_fastq_batch",
            &contig_names,
            &ordinals,
            &mate_offsets,
            &mate_indices,
            &reference_starts,
            &reference_ends,
            &sequences,
            &read_length,
            &quality,
            &paired_end
        )) {
        return NULL;
    }
    if (!PyTuple_Check(contig_names)) {
        PyErr_SetString(PyExc_ValueError, "FASTQ contig names must be a tuple");
        goto done;
    }
    if (read_length <= 0) {
        PyErr_SetString(PyExc_ValueError, "FASTQ read length must be positive");
        goto done;
    }
    if (quality < 33 || quality > 126) {
        PyErr_SetString(PyExc_ValueError, "FASTQ quality byte is outside Phred+33");
        goto done;
    }
    if (paired_end != 0 && paired_end != 1) {
        PyErr_SetString(PyExc_ValueError, "FASTQ paired-end flag must be boolean");
        goto done;
    }
    if (ordinals.len <= 0 || ordinals.len % 8 != 0) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ ordinals must contain packed uint64 values"
        );
        goto done;
    }
    fragment_count = ordinals.len / 8;
    if (PyTuple_GET_SIZE(contig_names) != fragment_count) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ contig names disagree with the fragment count"
        );
        goto done;
    }
    mates_per_fragment = paired_end ? 2 : 1;
    if (!checked_buffer_size(
            fragment_count,
            mates_per_fragment,
            &mate_count
        )) {
        goto done;
    }
    if (!checked_buffer_size(fragment_count + 1, 8, &expected_size)) {
        goto done;
    }
    if (mate_offsets.len != expected_size) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ mate offsets disagree with the fragment count"
        );
        goto done;
    }
    if (mate_indices.len != mate_count) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ mate indices disagree with the fragment count"
        );
        goto done;
    }
    if (!checked_buffer_size(mate_count, 8, &expected_size)) {
        goto done;
    }
    if (reference_starts.len != expected_size
        || reference_ends.len != expected_size) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ reference envelopes disagree with the mate count"
        );
        goto done;
    }
    if (!checked_buffer_size(mate_count, read_length, &expected_size)) {
        goto done;
    }
    if (sequences.len != expected_size) {
        PyErr_SetString(
            PyExc_ValueError,
            "FASTQ sequences disagree with the mate count and read length"
        );
        goto done;
    }
    ordinal_data = (const unsigned char *)ordinals.buf;
    offset_data = (const unsigned char *)mate_offsets.buf;
    mate_index_data = (const unsigned char *)mate_indices.buf;
    reference_start_data = (const unsigned char *)reference_starts.buf;
    reference_end_data = (const unsigned char *)reference_ends.buf;
    sequence_data = (const unsigned char *)sequences.buf;
    for (fragment_index = 0; fragment_index <= fragment_count; ++fragment_index) {
        const uint64_t expected_offset =
            (uint64_t)fragment_index * (uint64_t)mates_per_fragment;
        if (load_u64_le_bytes(offset_data + fragment_index * 8)
            != expected_offset) {
            PyErr_SetString(
                PyExc_ValueError,
                "FASTQ mate offsets are not canonical"
            );
            goto done;
        }
    }
    for (fragment_index = 0; fragment_index < fragment_count; ++fragment_index) {
        const Py_ssize_t mate_begin = fragment_index * mates_per_fragment;
        const uint64_t reference_start = load_u64_le_bytes(
            reference_start_data + mate_begin * 8
        );
        const uint64_t reference_end = load_u64_le_bytes(
            reference_end_data + mate_begin * 8
        );
        const char *contig_name;
        Py_ssize_t contig_name_length;
        Py_ssize_t mate_index;
        const int valid_name = valid_fastq_contig_name(
            PyTuple_GET_ITEM(contig_names, fragment_index),
            &contig_name,
            &contig_name_length
        );
        if (valid_name < 0) {
            goto done;
        }
        if (!valid_name) {
            PyErr_SetString(
                PyExc_ValueError,
                "FASTQ contig name is empty, non-printable, or contains whitespace"
            );
            goto done;
        }
        if (reference_start > reference_end || reference_start == UINT64_MAX) {
            PyErr_SetString(
                PyExc_ValueError,
                "FASTQ fragment reference envelope is invalid"
            );
            goto done;
        }
        for (mate_index = 1; mate_index < mates_per_fragment; ++mate_index) {
            const Py_ssize_t row = mate_begin + mate_index;
            if (load_u64_le_bytes(reference_start_data + row * 8)
                    != reference_start
                || load_u64_le_bytes(reference_end_data + row * 8)
                    != reference_end) {
                PyErr_SetString(
                    PyExc_ValueError,
                    "FASTQ mates disagree on the fragment reference envelope"
                );
                goto done;
            }
        }
    }
    for (fragment_index = 0; fragment_index < mate_count; ++fragment_index) {
        const unsigned char expected_index = (unsigned char)(
            fragment_index % mates_per_fragment
        );
        Py_ssize_t cycle;
        if (mate_index_data[fragment_index] != expected_index) {
            PyErr_SetString(
                PyExc_ValueError,
                "FASTQ mate rows are not ordered as 0 or 0,1"
            );
            goto done;
        }
        for (cycle = 0; cycle < read_length; ++cycle) {
            const unsigned char base =
                sequence_data[fragment_index * read_length + cycle];
            if (base != 'A' && base != 'C' && base != 'G'
                && base != 'T' && base != 'N') {
                PyErr_SetString(
                    PyExc_ValueError,
                    "FASTQ sequence contains a non-ACGTN byte"
                );
                goto done;
            }
        }
    }
    if (!checked_add_size(&record_base_size, 11, "FASTQ record")
        || read_length > (PY_SSIZE_T_MAX - record_base_size) / 2) {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_OverflowError, "FASTQ record exceeds Py_ssize_t");
        }
        goto done;
    }
    record_base_size += read_length * 2;
    record_lengths = PyTuple_New(fragment_count);
    if (record_lengths == NULL) {
        goto done;
    }
    for (fragment_index = 0; fragment_index < fragment_count; ++fragment_index) {
        const uint64_t ordinal = load_u64_le_bytes(
            ordinal_data + fragment_index * 8
        );
        const Py_ssize_t mate_begin = fragment_index * mates_per_fragment;
        const uint64_t reference_start = load_u64_le_bytes(
            reference_start_data + mate_begin * 8
        );
        const uint64_t reference_end = load_u64_le_bytes(
            reference_end_data + mate_begin * 8
        );
        const uint64_t left = reference_start + UINT64_C(1);
        const uint64_t right = reference_end > reference_start
            ? reference_end
            : left;
        const char *contig_name;
        Py_ssize_t contig_name_length;
        Py_ssize_t record_size = record_base_size;
        PyObject *lengths;
        contig_name = PyUnicode_AsUTF8AndSize(
            PyTuple_GET_ITEM(contig_names, fragment_index),
            &contig_name_length
        );
        if (contig_name == NULL) {
            goto done;
        }
        if (!checked_add_size(&record_size, contig_name_length, "FASTQ record")
            || !checked_add_size(
                &record_size,
                u64_decimal_length(left),
                "FASTQ record"
            )
            || !checked_add_size(
                &record_size,
                u64_decimal_length(right),
                "FASTQ record"
            )
            || !checked_add_size(
                &record_size,
                u64_decimal_length(ordinal),
                "FASTQ record"
            )) {
            goto done;
        }
        if (!checked_add_size(&read1_size, record_size, "FASTQ read1 batch")) {
            goto done;
        }
        if (paired_end
            && !checked_add_size(&read2_size, record_size, "FASTQ read2 batch")) {
            goto done;
        }
        lengths = Py_BuildValue(
            "(nnn)",
            record_size,
            paired_end ? record_size : 0,
            (Py_ssize_t)0
        );
        if (lengths == NULL) {
            goto done;
        }
        PyTuple_SET_ITEM(record_lengths, fragment_index, lengths);
    }
    read1 = PyBytes_FromStringAndSize(NULL, read1_size);
    if (read1 == NULL) {
        goto done;
    }
    if (paired_end) {
        read2 = PyBytes_FromStringAndSize(NULL, read2_size);
    } else {
        read2 = Py_None;
        Py_INCREF(read2);
    }
    if (read2 == NULL) {
        goto done;
    }
    read1_cursor = (unsigned char *)PyBytes_AS_STRING(read1);
    if (paired_end) {
        read2_cursor = (unsigned char *)PyBytes_AS_STRING(read2);
    }
    for (fragment_index = 0; fragment_index < fragment_count; ++fragment_index) {
        const uint64_t ordinal = load_u64_le_bytes(
            ordinal_data + fragment_index * 8
        );
        const Py_ssize_t mate_begin = fragment_index * mates_per_fragment;
        const uint64_t reference_start = load_u64_le_bytes(
            reference_start_data + mate_begin * 8
        );
        const uint64_t reference_end = load_u64_le_bytes(
            reference_end_data + mate_begin * 8
        );
        const char *contig_name;
        Py_ssize_t contig_name_length;
        contig_name = PyUnicode_AsUTF8AndSize(
            PyTuple_GET_ITEM(contig_names, fragment_index),
            &contig_name_length
        );
        if (contig_name == NULL) {
            goto done;
        }
        read1_cursor = write_fastq_record(
            read1_cursor,
            contig_name,
            contig_name_length,
            reference_start,
            reference_end,
            ordinal,
            1,
            sequence_data + mate_begin * read_length,
            read_length,
            (unsigned char)quality
        );
        if (paired_end) {
            read2_cursor = write_fastq_record(
                read2_cursor,
                contig_name,
                contig_name_length,
                reference_start,
                reference_end,
                ordinal,
                2,
                sequence_data + (mate_begin + 1) * read_length,
                read_length,
                (unsigned char)quality
            );
        }
    }
    if (
        read1_cursor
            != (unsigned char *)PyBytes_AS_STRING(read1) + read1_size
        || (paired_end
            && read2_cursor
                != (unsigned char *)PyBytes_AS_STRING(read2) + read2_size)
    ) {
        PyErr_SetString(PyExc_SystemError, "native FASTQ size accounting failed");
        goto done;
    }
    result = PyTuple_Pack(3, read1, read2, record_lengths);

done:
    Py_XDECREF(record_lengths);
    Py_XDECREF(read2);
    Py_XDECREF(read1);
    if (sequences.obj != NULL) {PyBuffer_Release(&sequences);}
    if (reference_ends.obj != NULL) {PyBuffer_Release(&reference_ends);}
    if (reference_starts.obj != NULL) {PyBuffer_Release(&reference_starts);}
    if (mate_indices.obj != NULL) {PyBuffer_Release(&mate_indices);}
    if (mate_offsets.obj != NULL) {PyBuffer_Release(&mate_offsets);}
    if (ordinals.obj != NULL) {PyBuffer_Release(&ordinals);}
    return result;
}


static int
initialize_site_attribute_names(void)
{
    if (column_site_index_name != NULL) {
        return 1;
    }
#define INTERN_COLUMN_NAME(target, value) \
    do { \
        target = PyUnicode_InternFromString(value); \
        if (target == NULL) { \
            return 0; \
        } \
    } while (0)
    INTERN_COLUMN_NAME(column_site_index_name, "site_index");
    INTERN_COLUMN_NAME(column_probability_name, "methylation_probability");
#undef INTERN_COLUMN_NAME
    return 1;
}


static PyObject *
required_cached_attribute(
    PyObject *record,
    PyObject *attribute_name,
    const char *display_name
)
{
    PyObject *value = PyObject_GetAttr(record, attribute_name);
    if (value == NULL) {
        PyErr_Format(
            PyExc_TypeError,
            "methylation site is missing required field %s",
            display_name
        );
    }
    return value;
}


static int
read_cached_unsigned_attribute(
    PyObject *record,
    PyObject *attribute_name,
    const char *display_name,
    uint64_t *result
)
{
    PyObject *value = required_cached_attribute(
        record,
        attribute_name,
        display_name
    );
    unsigned long long converted;
    if (value == NULL) {
        return 0;
    }
    converted = PyLong_AsUnsignedLongLong(value);
    Py_DECREF(value);
    if (converted == (unsigned long long)-1 && PyErr_Occurred()) {
        return 0;
    }
    *result = (uint64_t)converted;
    return 1;
}


static int
read_cached_double_attribute(
    PyObject *record,
    PyObject *attribute_name,
    const char *display_name,
    double *result
)
{
    PyObject *value = required_cached_attribute(
        record,
        attribute_name,
        display_name
    );
    double converted;
    if (value == NULL) {
        return 0;
    }
    converted = PyFloat_AsDouble(value);
    Py_DECREF(value);
    if (converted == -1.0 && PyErr_Occurred()) {
        return 0;
    }
    *result = converted;
    return 1;
}


static PyObject *
native_sample_bernoulli_sites(PyObject *self, PyObject *args)
{
    PyObject *sites;
    unsigned long long key;
    unsigned long long entity_ordinal;
    PyObject *result = NULL;
    Py_ssize_t index;
    Py_ssize_t count;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "O!KK:sample_bernoulli_sites",
            &PyTuple_Type,
            &sites,
            &key,
            &entity_ordinal
        )) {
        return NULL;
    }
    if (!initialize_site_attribute_names()) {
        return NULL;
    }
    count = PyTuple_GET_SIZE(sites);
    result = PyTuple_New(count);
    if (result == NULL) {
        return NULL;
    }
    for (index = 0; index < count; ++index) {
        PyObject *site = PyTuple_GET_ITEM(sites, index);
        uint64_t site_index;
        double probability;
        int sampled;
        PyObject *value;

        if (!read_cached_unsigned_attribute(
                site,
                column_site_index_name,
                "site_index",
                &site_index
            )
                || !read_cached_double_attribute(
                    site,
                    column_probability_name,
                    "methylation_probability",
                    &probability
                )) {
            Py_DECREF(result);
            return NULL;
        }
        if (!isfinite(probability) || probability < 0.0 || probability > 1.0) {
            Py_DECREF(result);
            PyErr_SetString(
                PyExc_ValueError,
                "methylation probability must be finite and in [0, 1]"
            );
            return NULL;
        }
        if (probability == 0.0) {
            sampled = 0;
        } else if (probability == 1.0) {
            sampled = 1;
        } else {
            const uint64_t random_value = philox_u64(
                (uint64_t)key,
                (uint64_t)entity_ordinal,
                site_index,
                0
            );
            const double uniform = ldexp((double)(random_value >> 11), -53);
            sampled = uniform < probability;
        }
        value = sampled ? Py_True : Py_False;
        Py_INCREF(value);
        PyTuple_SET_ITEM(result, index, value);
    }
    return result;
}


static PyMethodDef native_methods[] = {
    {
        "u64",
        native_u64,
        METH_VARARGS,
        PyDoc_STR("Return a Philox u64 for already validated arguments.")
    },
    {
        "bernoulli",
        native_bernoulli,
        METH_VARARGS,
        PyDoc_STR("Draw a Philox Bernoulli value for already validated arguments.")
    },
    {
        "philox_pairs",
        native_philox_pairs,
        METH_VARARGS,
        PyDoc_STR("Return packed Philox u64 pairs for validated counter buffers.")
    },
    {
        "format_fastq_batch",
        native_format_fastq_batch,
        METH_VARARGS,
        PyDoc_STR("Format validated columnar mates as exact FASTQ byte batches.")
    },
    {
        "sample_bernoulli_sites",
        native_sample_bernoulli_sites,
        METH_VARARGS,
        PyDoc_STR("Sample one validated fragment's Bernoulli methylation sites.")
    },
    {
        "apply_uniform_errors",
        native_apply_uniform_errors,
        METH_VARARGS,
        PyDoc_STR("Apply deterministic uniform sequencing errors to one mate.")
    },
    {
        "crc32c",
        native_crc32c,
        METH_VARARGS,
        PyDoc_STR("Return a Castagnoli CRC-32C for one contiguous bytes-like value.")
    },
    {
        "validate_protocol_batch_columns",
        bsreadsim_native_validate_protocol_batch_columns,
        METH_VARARGS,
        PyDoc_STR(
            "Semantically validate immutable protocol batch columns."
        )
    },
    {
        "pack_protocol_common_columns",
        bsreadsim_native_pack_protocol_common_columns,
        METH_VARARGS,
        PyDoc_STR(
            "Pack validated protocol common columns for NumPy processing."
        )
    },
    {
        "decode_protocol_fragments",
        bsreadsim_native_decode_protocol_fragments,
        METH_VARARGS,
        PyDoc_STR(
            "Reconstruct validated protocol Full-Truth typed fragments."
        )
    },
    {
        "canonical_truth_json_bytes",
        bsreadsim_native_truth_json_bytes,
        METH_VARARGS,
        PyDoc_STR(
            "Format one validated ProcessedFragment as canonical UTF-8 JSON bytes."
        )
    },
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef native_module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Private native hot loops for BSReadSim.",
    -1,
    native_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC
PyInit__native(void)
{
    initialize_crc32c_table();
    return PyModule_Create(&native_module);
}
