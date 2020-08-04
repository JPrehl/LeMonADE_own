## ###############  add your executables here and link them to LeMonADE ############# ##
add_executable(CrossLink1 lemonade_CrossLink_main.cpp)
target_link_libraries(CrossLink1 LeMonADE)

add_executable(CrossLinkSim lemonade_CrossLinkSim_main.cpp)
target_link_libraries(CrossLinkSim LeMonADE)
