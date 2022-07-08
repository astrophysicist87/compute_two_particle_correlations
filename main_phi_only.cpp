#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

//=================================================================
constexpr bool read_in_all_files = true;

constexpr size_t Dphi_bins = 32;
constexpr size_t n_mix = 10;

constexpr double pi = 3.1415926535897932384626433;
constexpr double twopi = 2.0*pi;
constexpr double Dphimin = -0.5*pi, Dphimax = 1.5*pi;
constexpr double Dphi_bw = (Dphimax-Dphimin)/Dphi_bins;

vector<double> phi_bin_centers, signal_pairs, mixed_pairs;

//=================================================================
size_t get_Dphi_bin( double Dphi_in )
{
	double Dphi = Dphi_in;
	while ( Dphi < Dphimin ) Dphi += twopi;
	while ( Dphi > Dphimax ) Dphi -= twopi;
	return static_cast<size_t>( (Dphi - Dphimin)/Dphi_bw );
}

//=================================================================
void load_arg_file(string filename, vector<string> & arguments);
void read_in_file( string filename, vector<double> & event );
void get_signal_pairs( const vector<double> & event );
void get_mixed_pairs( const vector<double> & event1, const vector<double> & event2 );

//=================================================================
int main(int argc, char *argv[])
{
  phi_bin_centers.resize(Dphi_bins);
	signal_pairs.resize(Dphi_bins);
	mixed_pairs.resize(Dphi_bins);

	for (size_t iDphibin = 0; iDphibin < Dphi_bins; iDphibin++)
    phi_bin_centers[iDphibin] = Dphimin + (iDphibin+0.5)*Dphi_bw;

	vector<string> arguments;
	if ( argc == 2 )
	{
		string argument_filename = argv[1];
		load_arg_file(argument_filename, arguments);
	}
	else
		for (size_t iArg = 1; iArg < argc; iArg++)
			arguments.push_back( argv[iArg] );

	const size_t nArguments = arguments.size();

//  for (string filename : arguments) cout << "Argument: " << filename << endl;
//  if (1) abort();


	vector<vector<double> > allEvents;
	if ( read_in_all_files )
//		for (size_t iArg = 0; iArg < nArguments; iArg++)
		for (string filename : arguments)
		{
cout << "Argument: " << filename << endl;
			vector<double> nextEvent;
//			string filename = arguments.at(iArg);
//      cout << "Currently using iArg = " << iArg << " in "
//            << arguments.size() << " vs " << argc << endl;
//			cout << "Reading in " << filename
//            << "(" << iArg+1 << " of " << nArguments << ")\n";
      cout << "Reading in " << filename << "\n";
			read_in_file( filename, nextEvent );
			allEvents.push_back( nextEvent );
		}

	// for choosing random events below
	vector<size_t> event_indices;
	for (size_t iArg = 0; iArg < nArguments; iArg++) event_indices.push_back( iArg );

	size_t count = 0;

	// loop over all files
	for (size_t iArg = 0; iArg < nArguments; iArg++)
	{
		string filename = arguments[iArg];

		const vector<double> & event = allEvents[iArg];

		// update signal pair distribution
		get_signal_pairs( event );

		// choose n_mix other random events to construct background
		std::vector<size_t> mix_events;
		std::sample(event_indices.begin(), event_indices.end(), 
                std::back_inserter(mix_events), n_mix + 1,	// extra event in case
                std::mt19937{std::random_device{}()});		// one is this event

		// loop over randomly chosen events and form background pairs
		size_t mixCount = 0;
		for ( const size_t & mix_event : mix_events )
		{
			if ( mixCount >= n_mix ) break;
			if ( mix_event == iArg ) continue;

			vector<double> & event_to_mix = allEvents[mix_event];

			get_mixed_pairs( event, event_to_mix );

			mixCount++;

		} // end loop over mixed events

		if (++count % 10 == 0) cout << "Completed " << count << " events so far!" << endl;

	} // end loop over all events

	// output ratio of signal to background
	cout << "Saving results to twoPC.dat" << endl;
	ofstream outfile("twoPC.dat");

	for (size_t iDphibin = 0; iDphibin < Dphi_bins; iDphibin++)
    outfile << phi_bin_centers[iDphibin] << "   "
            << signal_pairs[ iDphibin ] / ( mixed_pairs[ iDphibin ] + 1e-10 )
            << endl;

	outfile.close();

//double v_2_2 = signal_pairs[ iDphibin ] / ( mixed_pairs[ iDphibin ] + 1e-10 );

//cout << "v2{2} = " << v_2_2 << endl;

	return 0;
}

//==================================================================
void load_arg_file(string filename, vector<string> & arguments)
{
	string line;

	arguments.clear();
	ifstream infile( filename.c_str() );

	// read in arguments one line at a time
	if (infile.is_open())
		while ( getline (infile, line) )
			arguments.push_back( line );

	infile.close();
	return;
}

//==================================================================
void read_in_file(string filename, vector<double> & event)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, twopi);
	double random_phi = distribution(generator);

	event.clear();
	ifstream infile( filename.c_str() );
cout << "Reading in " << filename << endl;
	if (infile.is_open())
	{
		string line;
		double phi, eta; 
		while ( getline (infile, line) )
		{
			istringstream iss(line);
			iss >> phi >> eta;

			// impose cuts on eta
			if ( abs(eta) > 2.4 ) continue;
			event.push_back(phi + random_phi);
		}
cout << "event.size() = " << event.size() << endl;
	}
  else
  {
    cout << "Could not open " << filename << "!" << endl;
    abort();
  }

	infile.close();
	return;
}

//==================================================================
void get_signal_pairs( const vector<double> & event )
{
	const size_t event_size = event.size();

	for ( size_t i = 0;   i < event_size; i++ )
	for ( size_t j = i+1; j < event_size; j++ )
	{
    const size_t Dphi_ij = get_Dphi_bin( event[i]-event[j] );
    if ( Dphi_ij >= 0 && Dphi_ij < Dphi_bins )
      signal_pairs[ Dphi_ij ] += 1.0;
	}

	return;
}


//==================================================================
void get_mixed_pairs( const vector<double> & event1, const vector<double> & event2 )
{
	const size_t event1_size = event1.size();
	const size_t event2_size = event2.size();

	for ( size_t i = 0; i < event1_size; i++ )
	for ( size_t j = 0; j < event2_size; j++ )
	{
    const size_t Dphi_ij = get_Dphi_bin( event1[i]-event2[j] );
    mixed_pairs[ Dphi_ij ] += 1.0;
	}

	return;
}
