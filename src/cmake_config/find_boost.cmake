
# if(__find_maquis_boost)
#   return()
# endif()
# set(__find_maquis_boost YES)

function(find_maquis_boost)
	# parse arguments
	set(_nowhere)
	set(_curdest _nowhere)
	set(_val_args
		COMPONENTS)
	set(_bool_args
		REQUIRED)
	foreach(_arg ${_val_args} ${_bool_args})
		set(${_arg})
	endforeach()
	foreach(_element ${ARGN})
		list(FIND _val_args "${_element}" _val_arg_find)
		list(FIND _bool_args "${_element}" _bool_arg_find)
		if("${_val_arg_find}" GREATER "-1")
			set(_curdest "${_element}")
		elseif("${_bool_arg_find}" GREATER "-1")
			set("${_element}" ON)
			set(_curdest _nowhere)
		else()
			list(APPEND ${_curdest} "${_element}")
		endif()
	endforeach()

  message(STATUS "Looking for ${COMPONENTS}")
	if(_nowhere)
		message(FATAL_ERROR "Syntax error in use of find_maquis_boost!")
	endif()
  
  # exclude Boost Bindings from required
  set(FIND_Boost_REQUIRED)
  set(FIND_Boost_BINDINGS FALSE)
  set(FIND_Boost_TESTS FALSE)
  foreach(_package ${COMPONENTS})
    if(${_package} MATCHES "BINDINGS" OR ${_package} MATCHES "bindings")
      set(FIND_Boost_BINDINGS TRUE)
    else()
      if(${_package} MATCHES "unit_test")
        set(FIND_Boost_TESTS TRUE)
      endif(${_package} MATCHES "unit_test")
      list(APPEND FIND_Boost_REQUIRED ${_package})
    endif()
  endforeach()
  
  # Find libraries
  set(MAQUIS_Boost_FOUND MAQUIS_Boost_FOUND-NOTFOUND)
  set(MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY-NOTFOUND)
  if(ALPS_FOUND)
    set(MAQUIS_Boost_INCLUDE_DIRS ${ALPS_INCLUDE_DIRS})
    set(MAQUIS_Boost_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS})
    set(MAQUIS_Boost_LIBRARIES ${ALPS_Boost_LIBRARIES})
    
    if(FIND_Boost_TESTS)
      if(ALPS_HAS_BOOST_TEST)
    	  set(MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY "${ALPS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    	else(ALPS_HAS_BOOST_TEST)
    	  message(ERROR "ALPS not installing Boost.UnitTest")
    	endif(ALPS_HAS_BOOST_TEST)
  	endif(FIND_Boost_TESTS)
  	
  else(ALPS_FOUND)
    set(BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost sources)")
    find_package(Boost 1.43.0 COMPONENTS ${FIND_Boost_REQUIRED} REQUIRED)
    
    if(Boost_FOUND)
      set(MAQUIS_Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})  
      set(MAQUIS_Boost_LIBRARIES ${Boost_LIBRARIES})
      set(MAQUIS_Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
      if(FIND_Boost_TESTS)
        set(MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      endif(FIND_Boost_TESTS)
    endif(Boost_FOUND)
    
    if(FIND_Boost_BINDINGS)
      set(BOOST_BINDINGS_INCLUDE $ENV{BOOST_BINDINGS_INCLUDE} CACHE PATH "Path to Boost numeric bindings")
      list(APPEND MAQUIS_Boost_INCLUDE_DIRS ${BOOST_BINDINGS_INCLUDE})
    endif(FIND_Boost_BINDINGS)
  endif(ALPS_FOUND)
  
  set(BOOST_BINDINGS_INCLUDE "${BOOST_BINDINGS_INCLUDE}" CACHE PATH "Path to Boost numeric bindings" PARENT_SCOPE)
  set(MAQUIS_Boost_INCLUDE_DIRS "${MAQUIS_Boost_INCLUDE_DIRS}" PARENT_SCOPE)
  set(MAQUIS_Boost_LIBRARY_DIRS "${MAQUIS_Boost_LIBRARY_DIRS}" PARENT_SCOPE)
  set(MAQUIS_Boost_LIBRARIES "${MAQUIS_Boost_LIBRARIES}" PARENT_SCOPE)
  set(MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY "${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}" PARENT_SCOPE)
  
  
  
endfunction()
