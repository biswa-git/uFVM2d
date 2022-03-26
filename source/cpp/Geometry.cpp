#include<Geometry.hpp>

Geometry::Geometry() :m_is_read(false)
{

}

Geometry::~Geometry()
{
    /*
    if (!m_face_list.empty())
    {
        for (auto it : m_face_list)
        {
            FREE_OBJ_MACRO(*it);
        }
    }
    if (!m_vertex_list.empty())
    {
        for (auto it : m_vertex_list)
        {
            FREE_OBJ_MACRO(it);
        }
    }
    */
}

GeometryResult Geometry::Read(const std::string& file_name)
{
    GeometryResult result;

    if (!m_is_read)
    {
        std::ifstream file(file_name);
        if (file.is_open())
        {
            std::string line;
            std::istringstream iss;
            while (!file.eof())
            {
                getline(file, line);

                //-------------------------PHYSICAL NAMES---------------------------
                if (line == "$PhysicalNames")
                {
                    int number_of_physical_names, physical_name_type, physical_name_id;
                    std::string physical_name;

                    getline(file, line);
                    iss.clear();
                    iss.str(line);

                    if (!(iss >> number_of_physical_names)) { break; }

                    for (size_t i = 0; i < number_of_physical_names; i++)
                    {
                        getline(file, line);
                        iss.clear();
                        iss.str(line);
                        if (!(iss >> physical_name_type >> physical_name_id >> std::quoted(physical_name)))
                        {
                            break;
                        }

                        m_physical_group[physical_name_id] = PhysicalGroup::New(physical_name_type, physical_name);
                    }
                }

                //-------------------------NODES INFO---------------------------
                if (line == "$Nodes")
                {
                    int number_of_nodes, node_id;
                    double coords[3];

                    getline(file, line);
                    iss.clear();
                    iss.str(line);

                    if (!(iss >> number_of_nodes)) { break; }
                    for (size_t i = 0; i < number_of_nodes; i++)
                    {
                        getline(file, line);
                        iss.clear();
                        iss.str(line);

                        if (!(iss >> node_id >> coords[0] >> coords[1] >> coords[2]))
                        {
                            break;
                        }

                        m_vertex_list.emplace_back(Vertex::New(coords[0], coords[1], coords[2], node_id));
                    }
                }

                //-------------------------ELEMENTS INFO---------------------------
                if (line == "$Elements")
                {

                    long int num_of_elements, element_info[5], node_id[3];
                    Vector normal(0, 0, 1);
                    int offset = 0;
                    bool is_offset = false;
                    getline(file, line);
                    iss.clear();
                    iss.str(line);
                    if (!(iss >> num_of_elements)) { break; }
                    for (size_t i = 0; i < num_of_elements; i++)
                    {
                        getline(file, line);
                        iss.clear();
                        iss.str(line);


                        if (!(iss >> element_info[0] >> element_info[1] >> element_info[2] >> element_info[3] >> element_info[4]))
                        {
                            break;
                        }
                        
                        if (element_info[1] == PG_BOUNDARY)
                        {
                            if (!(iss >> node_id[0] >> node_id[1])) { break; };
                            m_physical_group[element_info[3]]->AddEdges(Edge::New(m_vertex_list[node_id[0] - 1], m_vertex_list[node_id[1] - 1]));
                        }
                        
                        if (element_info[1] == PG_INTERIOR)
                        {
                            if (!(iss >> node_id[0] >> node_id[1] >> node_id[2])) { break; }
                            if (!is_offset)
                            {
                                offset = element_info[0];
                                is_offset = true;
                            }
                            m_face_list.push_back
                            (
                                TriFace::New
                                (
                                    m_vertex_list[node_id[0] - 1],
                                    m_vertex_list[node_id[1] - 1],
                                    m_vertex_list[node_id[2] - 1],
                                    normal, element_info[0] - offset
                                )
                            );
                        }
                    }
                }
            }
        }

        //-------------------------INTERIOR EDGE---------------------------
        std::set<Edge*> interior_edge;
        for (auto face : m_face_list)
        {
            auto& half_edges = face->GetHalfEdge();

            for (auto half_edge : half_edges)
            {
                auto edge = half_edge->GetParentEdge();
                if (edge->GetHalfEdge(0)->GetFace() != nullptr && edge->GetHalfEdge(1)->GetFace() != nullptr)
                {
                    interior_edge.insert(edge);
                }
            }
        }

        for (auto edge: interior_edge)
        {
            m_interior_edge_list.push_back(edge);
        }

    }
    else
    {
        result.error = GEOMETRY_ALREADY_READ;
    }

    return result;
}

GeometryResult Geometry::Write(const std::string& file_location, const std::string& file_name)
{
    GeometryResult result;

    std::ofstream myfile;
    auto file = file_location + "/" + file_name + ".dat";
    myfile.open(file);
    myfile << "TITLE = \"title\"\n";
    myfile << "VARIABLES = \"X\", \"Y\", \"Z\"\n";
    myfile << "ZONE T = \"Rampant\", N = " << m_vertex_list.size() << ", E = " << m_face_list.size() << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n";

    for (auto it : m_vertex_list)
    {
        auto coord = it->GetPositionVector();
        myfile << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
    }

    for (auto it : m_face_list)
    {
        myfile << it->GetHalfEdge()[0]->GetStart()->GetId() << " " << it->GetHalfEdge()[1]->GetStart()->GetId() << " " << it->GetHalfEdge()[2]->GetStart()->GetId() << "\n";
    }
    
    myfile.close();

    result.success = true;
    return result;
}

std::vector<Vertex*>& Geometry::GetVertexList()
{
    return m_vertex_list;
}

const std::vector<Edge*>& Geometry::GetEdgeList() const
{
    return m_interior_edge_list;
}

const std::vector<Face*>& Geometry::GetFaceList() const
{
    return m_face_list;
}

std::map<int, PhysicalGroup*> Geometry::GetPhysicalGroup()
{
    return m_physical_group;
}

PhysicalGroup* Geometry::GetPhysicalgroupByName(std::string name)
{
    for (auto physical_group : m_physical_group)
    {
        if (physical_group.second->GetName() == name)
        {
            return physical_group.second;
        }
    }

    return nullptr;
}