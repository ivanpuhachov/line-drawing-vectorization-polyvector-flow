#include "FindAllPaths.h"
#include <stack>

std::vector<std::vector<vertex_descriptor>> findAllPaths(const G& g, vertex_descriptor s, vertex_descriptor t)
{
    std::vector<std::vector<vertex_descriptor>> result;
    std::stack<std::vector<vertex_descriptor>> partialPaths;
    partialPaths.push({s});
    
    while (!partialPaths.empty())
    {
        auto currentPath = partialPaths.top();
        partialPaths.pop();

        vertex_descriptor v = currentPath.back();

        for (auto [eit, e_end] = boost::out_edges(v, g); eit != e_end; ++eit)
        {
            if (eit->m_target == t)
            {
                auto newPath = currentPath;
                newPath.push_back(t);
                result.push_back(newPath);
            }
            else
            {
                if (std::find(currentPath.begin(), currentPath.end(), eit->m_target) == currentPath.end())
                {
                    auto newPath = currentPath;
                    newPath.push_back(eit->m_target);
                    partialPaths.push(newPath);
                }
            }
        }
    }

    std::cout << "Found " << result.size() << " paths: " << std::endl;
   
    for (int i=0; i<std::min((int)result.size(),20); ++i)
    {
        const auto& p = result[i];
        for (auto v : p)
        {
            std::cout << v << " ";
        }

        std::cout << std::endl;
    }
    return result;
}
