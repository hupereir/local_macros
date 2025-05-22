#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>

#include <regex>
#include <filesystem>

R__LOAD_LIBRARY(libodbc++.so)

// global connections
namespace
{
  std::unique_ptr<odbc::Connection> dbConnection;

  using file_list_t = std::vector<std::string>;
  std::ostream& operator << (std::ostream& out, const file_list_t& files )
  {
    for( const auto& file:files )
    { out << file << std::endl; }
    return out;
  }
}

//___________________________________________________________
file_list_t get_files( int runnumber )
{
  if( !dbConnection ) return file_list_t();

  // select filename from filelist  where runnumber = 53534 and hostname LIKE '%' || 'ebdc' || '%';
  const std::string sql = "SELECT filename FROM filelist WHERE runnumber = " + std::to_string(runnumber) + " and hostname LIKE '%' || 'ebdc' || '%';";
  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );

  file_list_t out;

  while (resultSet && resultSet->next())
  { out.emplace_back(resultSet->getString("filename")); }

  return out;
}


//___________________________________________________________
void GetRunFilesFromDB()
{
  try
  {
    dbConnection.reset( odbc::DriverManager::getConnection("daq", "", "") );
  }
  catch (odbc::SQLException &e)
  {
    std::cerr << "test_db_connection - Database connection failed: " << e.getMessage() << std::endl;
    return;
  }

  const int runnumber = 53534;
  auto files = get_files( runnumber );
  std::sort( files.begin(), files.end() );
  std::cout << "GetRunFilesFromDB - runnumber: " << runnumber << " files: " << files.size() << std::endl;

  file_list_t files_found;
  file_list_t files_missing;
  for( const auto& filename:files )
  {
    // destination path: /sphenix/lustre01/sphnxpro/physics/tpc/physics
    std::regex path_re("/bbox/bbox\\d/");
    const auto destination = std::regex_replace(filename, path_re, "/sphenix/lustre01/sphnxpro/physics/");

    if( std::filesystem::exists( destination ) )
    {
      std::cout << "Found " << destination << std::endl;
      files_found.emplace_back(destination);
    } else {
      files_missing.emplace_back(destination);
    }
  }

  std::sort(files_found.begin(), files_found.end());
  std::sort(files_missing.begin(), files_missing.end());

  std::cout << "GetRunFilesFromDB - runnumber: " << runnumber << " files_found: "
    << std::endl << files_found << std::endl;

  std::cout << "GetRunFilesFromDB - runnumber: " << runnumber << " files_missing: "
    << std::endl << files_missing << std::endl;

}
