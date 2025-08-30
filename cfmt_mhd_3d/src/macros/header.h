#define error_if         ERROR_IF_CPP(__FILE__, __LINE__)
#define raise_error(msg) RAISE_ERROR_CPP(__FILE__, __LINE__)(msg)

#define abort_if         ABORT_IF_CPP(__FILE__, __LINE__)
#define raise_abort(msg) RAISE_ABORT_CPP(__FILE__, __LINE__)(msg)

#define debug(msg)       DEBUG_CPP(__FILE__, __LINE__)(msg)

#define dump(x)          DUMP_CPP(__FILE__, __LINE__)(x)
