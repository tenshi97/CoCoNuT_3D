#include <stdio.h>
#include <stdlib.h>
/* Uncomment this to get rid of assert()
 * #define NDEBUG
 */
#include <assert.h>
#ifdef WITH_HDF5
#include <hdf5.h>
#endif

extern char _version_diff[];
extern size_t _version_diff_length;

extern char _version_log[];
extern size_t _version_log_length;

#ifdef HAVE_CONFIG
extern char _make_inc_config[];
extern size_t _make_inc_config_length;
#endif

extern char _version_commit_id[];
extern size_t _version_commit_id_length;

extern char _compile_time[];
extern size_t _compile_time_length;

#define _str(s) #s
#define str(s) _str(s)
const char *user = str(USER);
const char *machine = str(MACHINE);

void print_text_field(char *s, int num) {
        while (num--) putchar(*s++);
}

void version_short(void) {
	printf("Commit ");
	print_text_field(&_version_commit_id[0], _version_commit_id_length);
	printf("Compiled on %s by %s, ", machine, user);
	print_text_field(&_compile_time[0], _compile_time_length);
}

void version_log(void) {
        print_text_field(&_version_log[0], _version_log_length);
}

void version_diff(void) {
        print_text_field(&_version_diff[0], _version_diff_length);
}

#ifdef HAVE_CONFIG
void version_make_inc_config(void) {
	printf("# preprocessor configuration for this binary:\n");
        print_text_field(&_make_inc_config[0], _make_inc_config_length);
}
#endif

#ifdef WITH_HDF5
void c_write_string_dataset(hid_t loc_id, const char *name, const char *string, size_t len) {
	hid_t space_id, type_id, dset_id;

	/* Create scalar dataspace for the string attributes */
	space_id = H5Screate(H5S_SCALAR);
	assert(space_id >= 0);

	/* Create the datatype for this attribute */
	type_id = H5Tcopy(H5T_C_S1);
	assert(type_id >= 0);
	assert(H5Tset_size(type_id, len + 1) >= 0);

	/* Create the dataset */
	dset_id = H5Dcreate(loc_id, name, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//dset_id = H5Dcreate(loc_id, name, type_id, space_id, H5P_DEFAULT);
	assert(dset_id >= 0);

	/* Write the string */
	assert(H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, string) >= 0);

	/* Release resources */
	assert(H5Dclose(dset_id) >= 0);
	assert(H5Sclose(space_id) >= 0);
	assert(H5Tclose(type_id) >= 0);
}

void version_write_hdf5(void *locp) {
	hid_t loc = *((hid_t *) locp);
	c_write_string_dataset((hid_t) loc, "compile host", machine, sizeof(machine));
	c_write_string_dataset((hid_t) loc, "compile user", user, sizeof(user));
	c_write_string_dataset((hid_t) loc, "compile time", _compile_time, _compile_time_length);
	c_write_string_dataset((hid_t) loc, "commit id", _version_commit_id, _version_commit_id_length);
	c_write_string_dataset((hid_t) loc, "log entry", _version_log, _version_log_length);
	c_write_string_dataset((hid_t) loc, "diff to common base", _version_diff, _version_diff_length);
#ifdef HAVE_CONFIG
	c_write_string_dataset((hid_t) loc, "make.inc.vertex", _make_inc_config, _make_inc_config_length);
#endif
}
#endif
