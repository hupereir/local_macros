#include <RootUtil/FileManager.h>

#include <fun4all/Fun4AllUtils.h>

#include <fstream>
#include <iostream>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{

  // streamer for lists
  template <class T>
  std::ostream& operator<<(std::ostream& o, const std::vector<T>& list)
  {
    if (list.empty())
    {
      o << "{}";
    }
    else
    {

      int count = 0;

      const bool is_hex = (o.flags() & std::ios_base::hex);
      o << "  { " << std::endl << "    ";
      bool first = true;
      for (const auto& value : list)
      {
        if (count>0)
        {
          o << ", ";
        }

        if (is_hex)
        {
          o << "0x";
        }
        o << value;

        if( ++count == 10 )
        {
          count = 0;
          o << ", " << std::endl << "    ";
        }

        first = false;
      }
      o << std::endl << "  }";
    }
    return o;
  }


}

//____________________________________________
void CreateRunList()
{

  set_style(false);

  // get all processed runs and path
  std::map<int, TString> runmap;

  const std::vector<TString> path_list =
  {
    "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics",
    "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/fast",
    "/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/new_2024p002"
  };

  // create filemanager, load file selections
  for( const auto& path:path_list )
  {
    const auto selection = Form( "%s/run_*/*-00000.root", path.Data() );
    FileManager f( selection );

    // get single files
    const auto files = f.GetFiles();
    for( const auto& file:files )
    {
      const auto [runnumber,segment] = Fun4AllUtils::GetRunSegment(file.Data());
      if( runnumber )
      {
        runmap[runnumber]=path;
      } else {
        std::cout << "CheckRunHistory - could not get run number: " << file << std::endl;
      }
    }
  }

  std::cout << "CheckRunHistory - runmap: " << runmap.size() << std::endl;

  // write runs to file
  {

    // convert run map to list
    std::vector<int> runlist;
    runlist.reserve(runmap.size());
    for( const auto& [runnumber,path]:runmap ) { runlist.push_back( runnumber ); }

    const std::string runlist_filename( "macros/run_list_physics-new.h" );
    std::cout << "CreateRunList - runlist_filename: " << runlist_filename << std::endl;
    std::ofstream out( runlist_filename.c_str() );

    out << "#ifndef run_list_physics_h" << std::endl;
    out << "#define run_list_physics_h" << std::endl;
    out << "#include <vector>" << std::endl;
    out << "namespace RUNS" << std::endl;
    out << "{" << std::endl;
    out << "  std::vector<int> runnumbers = " << runlist << "; " << std::endl;
    out << "}" << std::endl;
    out << "#endif" << std::endl;
  }

  // generate jdf file
  {
    const std::string jdf_filename( "run_list_physics-new.jdf" );
    std::cout << "CreateRunList - jdf_filename: " << jdf_filename << std::endl;
    std::ofstream out( jdf_filename.c_str() );

    static constexpr int nevents = 500000;
    static constexpr int segment = 0;

    for( const auto& [runnumber, path]:runmap )
    {
      // out << "Arguments = " << path << " " << runnumber << " " << nevents << " " << segment << std::endl;
      out << "Arguments = " << nevents << " " << runnumber << " " << segment << " physics" << std::endl;
      out << "Queue 1" << std::endl;
      out << std::endl;
    }
  }
}
