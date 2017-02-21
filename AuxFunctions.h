#pragma once

#include <iostream>
#include <map>

namespace NuTo
{

    void WriteSimulationParameters(const std::map<std::string, std::string>& parameters, const std::string& fileName)
    {
        std::ofstream file(fileName + "parameters.txt");

        for (const auto& pair : parameters)
            file << pair.first << "\n" << pair.second << "\n\n";

        file.close();
    }

}
