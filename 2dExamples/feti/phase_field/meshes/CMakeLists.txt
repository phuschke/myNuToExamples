foreach(task RANGE 0 3)
    foreach(mesh wedge_split_test.mesh)

        create_symlink(${mesh}${task})

    endforeach()
endforeach()

foreach(task RANGE 0 4)
    foreach(mesh wedge_split_test_5_subdomains.mesh)

        create_symlink(${mesh}${task})

    endforeach()
endforeach()
