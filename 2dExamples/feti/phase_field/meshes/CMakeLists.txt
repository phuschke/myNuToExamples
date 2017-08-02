foreach(task RANGE 0 4)
    foreach(mesh wedge_split_test.mesh
            single_edge_notched_tension_test_quads.mesh
            single_edge_notched_tension_test_quads_coarse.mesh
            single_edge_notched_tension_test_quads_coarse_4_subs.mesh
            )
        create_symlink(${mesh}${task})
    endforeach()
endforeach()


foreach(task RANGE 0 4)
    foreach(mesh wedge_split_test_5_subdomains.mesh
            2d_single_edge_notched_tension_test_quads_5_subs.mesh)

        create_symlink(${mesh}${task})

    endforeach()
endforeach()

create_symlink(2d_single_edge_notched_tension_test_quads.msh)