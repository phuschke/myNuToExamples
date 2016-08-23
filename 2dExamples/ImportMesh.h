#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
using namespace std::string_literals;



struct Node
{
    Eigen::Vector3d mCoordinates;
    int mId;
};

struct Element
{
    int mId;
    std::vector<int> mNodeIds;
};

struct Boundary
{
    std::map<int, int> mNodeIdsMap;
    std::vector<int> mValues;
};

struct Interface
{
    std::map<int, int> mNodeIdsMap;
    int mValue;
};


struct ImportContainer
{
    using NodeList      = std::vector<Node>;
    using ElementList   = std::vector<Element>;
    using BoundaryList  = std::vector<Boundary>;
    using InterfaceList = std::vector<Interface>;

    NodeList        mNodeList;
    ElementList     mElementList;
    BoundaryList    mBoundaryList;
    InterfaceList   mInterfaceList;
};


std::vector<Node> ReadNodeData(std::ifstream& file)
{
    std::string line = "initialization";
    while(std::getline(file, line) and line.compare("Nodes")!=0);

    if (line.compare("Nodes")!=0)
        throw std::exception();

    int num_nodes = 0;
    file >> num_nodes;

    std::vector<Node> nodes(num_nodes);
    for (auto& node : nodes)
        file >> node.mId >> node.mCoordinates[0] >> node.mCoordinates[1] >> node.mCoordinates[2];

    return nodes;
}

std::vector<Element> ReadElementData(std::ifstream& file)
{
    std::string line = "initialization";
    while(std::getline(file, line) and line.compare("Elements")!=0);

    if (line.compare("Elements")!=0)
        throw std::exception();

    int num_elements = 0;
    file >> num_elements;

    std::vector<Element>   elements(num_elements);
    for (auto& element : elements)
    {
        file >> element.mId;

        element.mNodeIds.resize(4);

        for(auto& nodeId: element.mNodeIds)
            file >> nodeId;
    }

    return elements;
}

std::vector<Boundary> ReadBoundaryData(std::ifstream& file)
{
    std::string line = "initialization";
    while(std::getline(file, line) and line.compare("Boundaries")!=0);

    if (line.compare("Boundaries")!=0)
        throw std::exception();

    int num_boundaries = 0;
    file >> num_boundaries;

    std::vector<Boundary>   boundaries(num_boundaries);
    for (auto& boundary : boundaries)
    {
        int num_nodes = 0;
        file >> num_nodes;

        for(int i = 0; i < num_nodes; ++i)
        {
            int globalId  = 0;
            int localId   = 0;
            file >> globalId >> localId;
            boundary.mNodeIdsMap.emplace(globalId, localId);
        }

    }

    return boundaries;
}

std::vector<Interface> ReadInterfaceData(std::ifstream& file)
{
    std::string line = "initialization";
    while(std::getline(file, line) and line.compare("Interfaces")!=0);

    if (line.compare("Interfaces")!=0)
        throw std::exception();

    int num_interfaces = 0;
    file >> num_interfaces;

    std::vector<Interface>   interfaces(num_interfaces);
    for (auto& interface : interfaces)
    {
        int num_nodes = 0;
        file >> num_nodes;

        for(int i = 0; i < num_nodes; ++i)
        {
            int globalId  = 0;
            int localId   = 0;
            file >> globalId >> localId;
            interface.mNodeIdsMap.emplace(globalId, localId);
        }

        file >> interface.mValue;


    }

    return interfaces;
}


ImportContainer ImportMeshFile(const std::string &rFileName)
{

    std::ifstream file(rFileName.c_str(), std::ios::in);
    if (not file.is_open())
        throw std::exception();

    ImportContainer importContainer;

    importContainer.mNodeList         = ReadNodeData              (file);
    importContainer.mElementList      = ReadElementData           (file);
    importContainer.mBoundaryList     = ReadBoundaryData          (file);
    importContainer.mInterfaceList    = ReadInterfaceData         (file);


    file.close();

    return importContainer;
}
