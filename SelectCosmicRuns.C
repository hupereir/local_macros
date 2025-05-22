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

  struct run_info_t
  {
    int runnumber = 0;
    uint64_t events = 0;
    bool has_tpot = false;
    bool has_tpc = false;
    bool has_magnet = false;
    bool selected() const { return has_tpot && has_tpc; }
  };

  using run_info_list_t = std::vector<run_info_t>;

  std::ostream& operator << (std::ostream& out, const run_info_t& info)
  {
    out
      << "runnumber: " << info.runnumber
      << " events: " << info.events
      << " has_tpot: " << info.has_tpot
      << " has_tpc: " << info.has_tpc
      << " has_magnet: " << info.has_magnet;

    return out;
  }

  std::ostream& operator << (std::ostream& out, const run_info_list_t& list)
  {
    for( const auto& info:list)
    {
      out << info << std::endl;
    }
    return out;
  }
}

//_________________________________________
bool connect()
{

 try
  {
    dbConnection.reset( odbc::DriverManager::getConnection("daq", "", "") );
  }
  catch (odbc::SQLException &e)
  {
    std::cerr << "test_db_connection - Database connection failed: " << e.getMessage() << std::endl;
    return false;
  }
  return true;
}

//_________________________________________
bool has_tpot( const int runnumber )
{
  const std::string sql =
    "SELECT filename FROM filelist "
    "WHERE runnumber = " + std::to_string(runnumber) +
    "AND hostname LIKE '%' || 'ebdc39' || '%';";

  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );
  return resultSet && resultSet->next();
}

//_________________________________________
bool has_tpc( const int runnumber )
{
  const std::string sql =
    "SELECT filename FROM filelist "
    "WHERE runnumber = " + std::to_string(runnumber) +
    "AND hostname LIKE '%' || 'ebdc02' || '%';";

  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );
  return resultSet && resultSet->next();
}


//_________________________________________
bool has_magnet( const int runnumber )
{
  const std::string sql =
    "SELECT magnet_on FROM magnet_info "
    "WHERE runnumber = " + std::to_string(runnumber);
  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );
  return resultSet && resultSet->next() && resultSet->getBoolean("magnet_on");
}

//_________________________________________
run_info_list_t selectRuns()
{
  const std::string sql =
    "SELECT runnumber, eventsinrun FROM run "
    "WHERE runtype='cosmics' "
    "AND runnumber>0 "
    "AND EXTRACT(EPOCH FROM (ertimestamp - brtimestamp)) > 1800 "
    "AND brtimestamp>'2024-01-01'";

  std::unique_ptr<odbc::Statement> stmt( dbConnection->createStatement() );
  std::unique_ptr<odbc::ResultSet> resultSet( stmt->executeQuery(sql) );

  run_info_list_t run_info_list;
  while( resultSet && resultSet->next() )
  {
    run_info_t info;
    info.runnumber = resultSet->getInt("runnumber");
    info.events = resultSet->getInt("eventsinrun");
    run_info_list.emplace_back(info);
  }

  return run_info_list;
}

void SelectCosmicRuns()
{
  connect();
  auto info_list = selectRuns();
  for( auto&& info:info_list )
  {
    info.has_tpot = has_tpot(info.runnumber);
    info.has_tpc = has_tpc(info.runnumber);
    info.has_magnet = has_magnet(info.runnumber);
  }

  std::cout << info_list << std::endl;


  static constexpr float occupancy = 2e-4;

  {
    // total runs with TPOT
    auto n_runnumbers = std::count_if( info_list.begin(), info_list.end(),
      []( const run_info_t& info )
    { return info.selected(); } );

    auto n_events = std::accumulate(  info_list.begin(), info_list.end(), (uint64_t) 0,
      []( uint64_t result, const run_info_t info ) { return info.selected() ? result+info.events:result; } );

    std::cout << "Selected run count: " << n_runnumbers << std::endl;
    std::cout << "Selected event count: " << Form( "%.3g", float(n_events)) << std::endl;
    std::cout << "Selected track counts = " << Form( "%.3g", occupancy*n_events) << std::endl;
  }

  // magnets on runs
  {
    // total runs with TPOT
    auto n_runnumbers = std::count_if( info_list.begin(), info_list.end(),
      []( const run_info_t& info )
    { return info.selected()&&info.has_magnet; } );

    auto n_events = std::accumulate(  info_list.begin(), info_list.end(), (uint64_t) 0,
      []( uint64_t result, const run_info_t info ) { return info.selected()&&info.has_magnet ? result+info.events:result; } );

    std::cout << "Selected magnet on run count: " << n_runnumbers << std::endl;
    std::cout << "Selected magnet on event count: " << Form( "%.3g", float(n_events)) << std::endl;
    std::cout << "Selected magnet on track counts = " << Form( "%.3g", occupancy*n_events) << std::endl;
  }

}
