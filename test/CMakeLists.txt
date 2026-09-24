SET(TEST_WORK_ROOT ${CMAKE_CURRENT_BINARY_DIR}/work)
SET(TEST_DATA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_LIBRARY(REGRESSION_UTILS OBJECT src/REGRESSION_UTILS.f90)

ADD_EXECUTABLE(dmrg_test_spin
  src/Spin1d_SU2/spinSU2_1d.f90
  $<TARGET_OBJECTS:REGRESSION_UTILS>
)
ADD_EXECUTABLE(dmrg_test_hubbard
  src/Hubbard1d/hm_1d.f90
  $<TARGET_OBJECTS:REGRESSION_UTILS>
)

ADD_DEPENDENCIES(dmrg_test_spin dmrg REGRESSION_UTILS)
ADD_DEPENDENCIES(dmrg_test_hubbard dmrg REGRESSION_UTILS)
TARGET_LINK_LIBRARIES(dmrg_test_spin dmrg)
TARGET_LINK_LIBRARIES(dmrg_test_hubbard dmrg)

FOREACH(model Spin1d_SU2 Hubbard1d)
  SET(work_dir ${TEST_WORK_ROOT}/${model})
  FILE(REMOVE_RECURSE ${work_dir})
  FILE(MAKE_DIRECTORY ${work_dir})
  IF(model STREQUAL "Spin1d_SU2")
    SET(check_files energy.check entropy.check spin_local.check
        spin_nn.check spin_1j.check)
    CONFIGURE_FILE(${TEST_DATA_ROOT}/${model}/DMRG.conf.SpinSU2
                   ${work_dir}/DMRG.conf COPYONLY)
  ELSEIF(model STREQUAL "Hubbard1d")
    SET(check_files energy.check entropy.check hubbard_local.check
        density_nn.check density_1j.check)
    CONFIGURE_FILE(${TEST_DATA_ROOT}/${model}/DMRG.conf.Hubbard
                   ${work_dir}/DMRG.conf COPYONLY)
  ENDIF()
  #Stage only the approved reference tables required by each test driver.
  #The files in test/reference are used only to generate these data manually.
  FOREACH(check_file IN LISTS check_files)
    CONFIGURE_FILE(${TEST_DATA_ROOT}/${model}/${check_file}
                   ${work_dir}/${check_file} COPYONLY)
  ENDFOREACH()
ENDFOREACH()

FUNCTION(add_dmrg_test test_name target model)
  SET(work_dir ${TEST_WORK_ROOT}/${model})
  IF(USE_MPI)
    ADD_TEST(NAME ${test_name}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target}>
              FINPUT=DMRG.conf ${MPIEXEC_POSTFLAGS})
    SET_TESTS_PROPERTIES(${test_name} PROPERTIES PROCESSORS 2)
  ELSE()
    ADD_TEST(NAME ${test_name}
      COMMAND $<TARGET_FILE:${target}> FINPUT=DMRG.conf)
  ENDIF()
  SET_TESTS_PROPERTIES(${test_name} PROPERTIES
    WORKING_DIRECTORY ${work_dir}
    LABELS "regression;dmrg;${model}"
    TIMEOUT 1200
    ENVIRONMENT "OMP_NUM_THREADS=1"
  )
ENDFUNCTION()

add_dmrg_test(dmrg_spin_regression dmrg_test_spin Spin1d_SU2)
add_dmrg_test(dmrg_hubbard_regression dmrg_test_hubbard Hubbard1d)
