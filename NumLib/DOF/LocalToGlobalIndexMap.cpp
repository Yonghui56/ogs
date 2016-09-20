/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalToGlobalIndexMap.h"

namespace NumLib
{

template <typename ElementIterator>
void
LocalToGlobalIndexMap::findGlobalIndices(
    ElementIterator first, ElementIterator last,
    std::size_t const mesh_id,
    const unsigned comp_id, const unsigned comp_id_write)
{
    _rows.resize(std::distance(first, last), _mesh_subsets.size());

    // For each element find the global indices for node/element
    // components.
    std::size_t elem_id = 0;
    for (ElementIterator e = first; e != last; ++e, ++elem_id)
    {
        std::size_t const nnodes = (*e)->getNumberOfNodes();

        LineIndex indices;
        indices.reserve(nnodes);

        for (unsigned n = 0; n < nnodes; n++)
        {
            MeshLib::Location l(mesh_id,
                                MeshLib::MeshItemType::Node,
                                (*e)->getNode(n)->getID());
            indices.push_back(_mesh_component_map.getGlobalIndex(l, comp_id));
        }

        _rows(elem_id, comp_id_write) = std::move(indices);
    }
}


LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>>&& mesh_subsets,
    NumLib::ComponentOrder const order)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(_mesh_subsets, order)
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.


    _variable_component_offsets.reserve(_mesh_subsets.size());
    std::size_t offset = 0;
    for (int variable_id = 0; variable_id < _mesh_subsets.size();
         ++variable_id)
    {
        _variable_component_offsets.push_back(offset);
        auto const& mss = *_mesh_subsets[variable_id];
        for (int component_id = 0; component_id < mss.size();
             ++component_id)
        {
            MeshLib::MeshSubset const& ms = mss.getMeshSubset(component_id);
            std::size_t const mesh_id = ms.getMeshID();

            auto const global_component_id =
                getGlobalComponent(variable_id, component_id);

            findGlobalIndices(ms.elementsBegin(), ms.elementsEnd(), mesh_id,
                              global_component_id, global_component_id);
        }

        // increase by number of components of that variable
        offset += mss.size();
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>>&& mesh_subsets,
    int const component_id,
    std::vector<MeshLib::Element*> const& elements,
    NumLib::MeshComponentMap&& mesh_component_map)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(std::move(mesh_component_map)),
      _variable_component_offsets(1, 0) // Single variable only.
{
    // There is only on mesh_subsets in the vector _mesh_subsets.
    assert(_mesh_subsets.size() == 1);
    auto const mss = *_mesh_subsets.front();

    // For all MeshSubset in mesh_subsets and each element of that MeshSubset
    // save a line of global indices.
    for (MeshLib::MeshSubset const* const ms : mss)
    {
        std::size_t const mesh_id = ms->getMeshID();

        findGlobalIndices(elements.cbegin(), elements.cend(), mesh_id,
                          component_id, 0);  // There is only one component to
                                             // write out, therefore the zero
                                             // parameter.
    }
}

LocalToGlobalIndexMap* LocalToGlobalIndexMap::deriveBoundaryConstrainedMap(
    int const variable_id,
    int const component_id,
    std::unique_ptr<MeshLib::MeshSubsets>&& mesh_subsets,
    std::vector<MeshLib::Element*> const& elements) const
{
    DBUG("Construct reduced local to global index map.");
    // Create a subset of the current mesh component map.
    auto const global_component_id =
        getGlobalComponent(variable_id, component_id);

    auto mesh_component_map =
        _mesh_component_map.getSubset(global_component_id, *mesh_subsets);

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    all_mesh_subsets.emplace_back(std::move(mesh_subsets));
    return new LocalToGlobalIndexMap(std::move(all_mesh_subsets),
                                     global_component_id, elements,
                                     std::move(mesh_component_map));
}

std::size_t
LocalToGlobalIndexMap::dofSizeWithGhosts() const
{
    return _mesh_component_map.dofSizeWithGhosts();
}

std::size_t
LocalToGlobalIndexMap::size() const
{
    return _rows.rows();
}

LocalToGlobalIndexMap::RowColumnIndices
LocalToGlobalIndexMap::operator()(std::size_t const mesh_item_id, const unsigned component_id) const
{
    return RowColumnIndices(_rows(mesh_item_id, component_id),
                            _columns(mesh_item_id, component_id));
}

std::size_t
LocalToGlobalIndexMap::getNumElementDOF(std::size_t const mesh_item_id) const
{
    std::size_t ndof = 0;

    for (unsigned c=0; c<_rows.cols(); ++c)
    {
        ndof += _rows(mesh_item_id, c).size();
    }

    return ndof;
}

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map)
{
    std::size_t const max_lines = 10;
    std::size_t lines_printed = 0;

    os << "Rows of the local to global index map; " << map._rows.size()
        << " rows\n";
    for (std::size_t e=0; e<map.size(); ++e)
    {
        for (std::size_t c=0; c<map.getNumComponents(); ++c)
        {
            auto const& line = map._rows(e, c);

            os << "c" << c << " { ";
            std::copy(line.cbegin(), line.cend(),
                std::ostream_iterator<std::size_t>(os, " "));
            os << " }\n";
        }

        if (lines_printed++ > max_lines)
        {
            os << "...\n";
            break;
        }
    }
    lines_printed = 0;

    os << "Mesh component map:\n" << map._mesh_component_map;
    return os;
}
#endif  // NDEBUG

}   // namespace NumLib
