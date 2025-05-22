#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>

R__LOAD_LIBRARY(libodbc++.so)

#include "run_list_physics.h"

// global connections
namespace
{
  std::unique_ptr<odbc::Connection> dbConnection;
}

//___________________________________________________________
int get_events( int runnumber )
{
  if( !dbConnection ) return 0;

  // select nr_events from event_numbers where runnumber = 52078 and hostname='gl1daq'
  const std::string sql = "SELECT nr_events FROM event_numbers WHERE runnumber = " + std::to_string(runnumber) + " and hostname='gl1daq';";
  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );

  if (resultSet && resultSet->next())
  { return resultSet->getInt("nr_events"); }

  return 0;
}


//___________________________________________________________
void GetRunEventsFromDB()
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

  for( const auto& runnumber:runnumbers )
  {
    std::cout << "GetRunEventsFromDB - runnumber: " << runnumber << " events: " << get_events(runnumber) << std::endl;
  }

}
