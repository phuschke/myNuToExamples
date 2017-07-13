
#include <jsoncpp/json/json.h>
#include <iostream>
#include <fstream>

int main()
{
    Json::Value root;
    Json::Reader reader;


    std::ifstream file("mesh.json");

    if (not reader.parse(file, root, false))
        std::cout << "Oh no!" << std::endl;


    std::cout << "Hello " << std::endl;


    std::vector<int> vec;

    for (const auto& bla : root["LocalToGlobalMap"])
        vec.push_back(bla.asInt());

    for (const auto& bla : vec)
        std::cout << bla << std::endl;


    std::cout << root["Elements"][0]["NodalConnectivity"][0][0] << std::endl;


    std::cout << root["Elements"][0]["NodalConnectivity"].size() << std::endl;
}
