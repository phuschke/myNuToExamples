#pragma once

#include <iostream>
#include <map>

namespace NuTo
{

    void WriteSimulationParameters(const std::map<std::string, std::string>& parameters, const std::string& fileName)
    {
        std::cout << "Writing parameters to file: " << fileName + "parameters.txt" << std::endl;
        std::ofstream file(fileName + "parameters.txt");

        if (not file.is_open())
        {
            std::cout << "cant open file: " << fileName << std::endl;
        }

        for (const auto& pair : parameters)
            file << pair.first << "\n" << pair.second << "\n\n";

        file.close();
    }

}
