# Writes OUT (a header defining MEOSAR_VERSION), only when the version changed.
# The version: the latest v* tag as git describes it (v0.1.0, or v0.1.0-3-gabc1234 three commits
# on), FALLBACK+gHASH before the first tag, FALLBACK outside a git checkout; -dirty with changes.
execute_process(COMMAND git -C "${SRC}" describe --tags --match "v[0-9]*" --dirty
                OUTPUT_VARIABLE V OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE rc)
if(rc EQUAL 0 AND V MATCHES "^v")
  string(SUBSTRING "${V}" 1 -1 V)
else()
  execute_process(COMMAND git -C "${SRC}" describe --always --dirty
                  OUTPUT_VARIABLE H OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE rc)
  if(rc EQUAL 0 AND NOT H STREQUAL "")
    set(V "${FALLBACK}+g${H}")
  else()
    set(V "${FALLBACK}")
  endif()
endif()
set(TXT "#define MEOSAR_VERSION \"${V}\"\n")
set(OLD "")
if(EXISTS "${OUT}")
  file(READ "${OUT}" OLD)
endif()
if(NOT OLD STREQUAL TXT)
  file(WRITE "${OUT}" "${TXT}")
endif()
