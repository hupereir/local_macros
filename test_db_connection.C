#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>

R__LOAD_LIBRARY(libodbc++.so)

// based on code from

void test_db_connection( int runnumber = 52078)
{
  std::unique_ptr<odbc::Connection> dbConnection;
  try
  {
    dbConnection.reset( odbc::DriverManager::getConnection("daq", "", "") );
  }
  catch (odbc::SQLException &e)
  {
    std::cerr << "test_db_connection - Database connection failed: " << e.getMessage() << std::endl;
    return;
  }

  // select nr_events from event_numbers where runnumber = 52078 and hostname='gl1daq'
  const std::string sql = "SELECT nr_events FROM event_numbers WHERE runnumber = " + std::to_string(runnumber) + " and hostname='gl1daq';";
  odbc::Statement *stmt = dbConnection->createStatement();
  odbc::ResultSet *resultSet = stmt->executeQuery(sql);

  if (resultSet && resultSet->next())
  {
    const auto events = resultSet->getInt("nr_events");
    std::cout << "test_db_connection - runnumber: " << runnumber << " events: " << events << std::endl;
  }

}
