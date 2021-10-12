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
constexpr size_t Deta_bins = 28;
constexpr size_t n_mix = 10;

constexpr double pi = 3.1415926535897932384626433;
constexpr double twopi = 2.0*pi;
constexpr double Dphimin = -0.5*pi, Dphimax = 1.5*pi;
constexpr double Detamin = -4.0, Detamax = 4.0;
constexpr double Dphi_bw = (Dphimax-Dphimin)/Dphi_bins;
constexpr double Deta_bw = (Detamax-Detamin)/Deta_bins;

vector<double> signal_pairs, mixed_pairs;

//=================================================================
inline size_t indexer( size_t i, size_t j )
{
	return i * Deta_bins + j;
}

//=================================================================
size_t get_Dphi_bin( double Dphi_in )
{
	double Dphi = Dphi_in;
	while ( Dphi < Dphimin ) Dphi += twopi;
	while ( Dphi > Dphimax ) Dphi -= twopi;
	return static_cast<size_t>( (Dphi - Dphimin)/Dphi_bw );
}

//=================================================================
size_t get_Deta_bin( double Deta )
{
	return static_cast<size_t>( (Deta - Detamin)/Deta_bw );
}

//=================================================================
void load_arg_file(string filename, vector<string> & arguments);
void read_in_file( string filename, vector<double> & event );
void get_signal_pairs( const vector<double> & event );
void get_mixed_pairs( const vector<double> & event1, const vector<double> & event2 );

//=================================================================
int main(int argc, char *argv[])
{
	signal_pairs.resize(Dphi_bins*Deta_bins);
	mixed_pairs.resize(Dphi_bins*Deta_bins);

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

	vector<vector<double> > allEvents;
	if ( read_in_all_files )
		for (size_t iArg = 0; iArg < nArguments; iArg++)
		{
			vector<double> nextEvent;
			string filename = arguments[iArg];
			cout << "Reading in " << filename << "\n";
			read_in_file( filename, nextEvent );
			allEvents.push_back( nextEvent );
		}

	// for choosing random events below
	vector<size_t> event_indices;
	for (size_t iArg = 0; iArg < nArguments; iArg++) event_indices.push_back( iArg );

//	// try generating all mix events at once
//	for (size_t iArg = 0; iArg < nArguments; iArg++)
//	{
//		// choose n_mix other random events to construct background
//		std::vector<size_t> mix_events;
//		std::sample(event_indices.begin(), event_indices.end(), 
//				std::back_inserter(mix_events), n_mix + 1,	// extra event in case
//				std::mt19937{std::random_device{}()});		// one is this event
//		all_mix_events.push_back( mix_events );
//		cout << "Generated mix_events for iArg = " << iArg << " \n";
//	}




	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////

	vector<vector<size_t> > all_mix_events;
	vector<vector<vector<size_t> > > all_mix_events_per_thread( omp_get_max_threads() );
	const size_t arguments_per_thread = nArguments / omp_get_max_threads();
	#pragma omp parallel for
	for (size_t iThread = 0; iThread < omp_get_max_threads(); iThread++)
	{
		vector<size_t> event_indices_local = event_indices;
		vector<vector<size_t> > & all_mix_events_local = all_mix_events_per_thread[iThread];
		all_mix_events_local.resize( arguments_per_thread );

		// try generating all mix events at once
		for (size_t iArg = iThread * arguments_per_thread;
					iArg < (iThread+1) * arguments_per_thread; iArg++)
		{
			// choose n_mix other random events to construct background
			std::vector<size_t> mix_events;
			std::sample(event_indices_local.begin(), event_indices_local.end(), 
					std::back_inserter(mix_events), n_mix + 1,	// extra event in case
					std::mt19937{std::random_device{}()});		// one is this event
			all_mix_events_local[ iArg - iThread * arguments_per_thread ] = mix_events;
			cout << "check: " << iThread << "   " << iArg << "\n";
		}
	}

	for (size_t iThread = 0; iThread < omp_get_max_threads(); iThread++)
	for (const auto & events_to_mix_here : all_mix_events_per_thread[iThread] )
		all_mix_events.push_back( events_to_mix_here );


	while ( all_mix_events.size() < nArguments )
	{
		// choose n_mix other random events to construct background
		std::vector<size_t> mix_events;
		std::sample(event_indices.begin(), event_indices.end(), 
				std::back_inserter(mix_events), n_mix + 1,	// extra event in case
				std::mt19937{std::random_device{}()});		// one is this event
		all_mix_events.push_back( mix_events );
	}

	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////


	size_t count = 0;

	// loop over all files
	for (size_t iArg = 0; iArg < nArguments; iArg++)
	{
		//cout << "Reading in " << arguments[iArg] << "; ";

		string filename = arguments[iArg];
		/*vector<double> event;
		if (!read_in_all_files)
			read_in_file( filename, event );
		else
			event = allEvents[iArg];*/
		const vector<double> & event = allEvents[iArg];

		//cout << "read in " << event.size()/2 << " particles.\n";

		// update signal pair distribution
		//cout << "\t - updating signal distribution...";
		get_signal_pairs( event );
		//cout << "done!" << endl;

		/*
		// choose n_mix other random events to construct background
		std::vector<size_t> mix_events;
		std::sample(event_indices.begin(), event_indices.end(), 
                std::back_inserter(mix_events), n_mix + 1,	// extra event in case
                std::mt19937{std::random_device{}()});		// one is this event
		*/
		const std::vector<size_t> & mix_events = all_mix_events[iArg];

		// loop over randomly chosen events and form background pairs
		size_t mixCount = 0;
		for ( const size_t & mix_event : mix_events )
		{
			if ( mixCount >= n_mix ) break;
			if ( mix_event == iArg ) continue;

			//cout << "\t - mixing " << arguments[iArg] << " with " << arguments[mix_event] << "\n";

			/*vector<double> event_to_mix;
			if (!read_in_all_files)
				read_in_file( arguments[mix_event], event_to_mix );
			else
				event_to_mix = allEvents[mix_event];*/
			vector<double> & event_to_mix = allEvents[mix_event];

			//cout << "\t - updating mixed distribution...";
			get_mixed_pairs( event, event_to_mix );
			//cout << "done!" << "\n";

			mixCount++;

		} // end loop over mixed events

		if (++count % 1000 == 0) cout << "Completed " << count << " events so far!" << endl;

	} // end loop over all events

	//cout << flush;

	// output ratio of signal to background
	cout << "Saving results to twoPC.dat" << endl;
	ofstream outfile("twoPC.dat");
	for (size_t iDphibin = 0; iDphibin < Dphi_bins; iDphibin++)
	{
		for (size_t iDetabin = 0; iDetabin < Deta_bins; iDetabin++)
			outfile << signal_pairs[ indexer( iDphibin, iDetabin ) ]
						/ ( mixed_pairs[ indexer( iDphibin, iDetabin ) ] + 1e-10 ) << "   ";
		outfile << endl;
	}
	//outfile << endl;
	outfile.close();

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
			event.push_back(eta);
		}
	}

	infile.close();
	return;
}

//==================================================================
void get_signal_pairs( const vector<double> & event )
{
	const size_t event_size = event.size()/2;

	for ( size_t i = 0;   i < event_size; i++ )
	for ( size_t j = i+1; j < event_size; j++ )
	{
		const size_t Deta_ij = get_Deta_bin(event[2*i+1]-event[2*j+1]);
		if ( Deta_ij >= 0 && Deta_ij < Deta_bins )
		{
			const size_t Dphi_ij = get_Dphi_bin(event[2*i+0]-event[2*j+0]);
			if ( Dphi_ij >= 0 && Dphi_ij < Dphi_bins )
				signal_pairs[ indexer( Dphi_ij, Deta_ij ) ] += 1.0;
		}

		// ensure reflection symmetry under eta --> -eta (is this correct?)
		const size_t Deta_ji = get_Deta_bin(event[2*j+1]-event[2*i+1]);
		if ( Deta_ji >= 0 && Deta_ji < Deta_bins )
		{
			const size_t Dphi_ji = get_Dphi_bin(event[2*j+0]-event[2*i+0]);
			if ( Dphi_ji >= 0 && Dphi_ji < Dphi_bins )
				signal_pairs[ indexer( Dphi_ji, Deta_ji ) ] += 1.0;
		}
	}

	return;
}


//==================================================================
void get_mixed_pairs( const vector<double> & event1, const vector<double> & event2 )
{
	const size_t event1_size = event1.size()/2;
	const size_t event2_size = event2.size()/2;

	for ( size_t i = 0; i < event1_size; i++ )
	for ( size_t j = 0; j < event2_size; j++ )
	{
		const size_t Deta_ij = get_Deta_bin(event1[2*i+1]-event2[2*j+1]);
		if ( Deta_ij >= 0 && Deta_ij < Deta_bins )
		{
			const size_t Dphi_ij = get_Dphi_bin(event1[2*i+0]-event2[2*j+0]);
			mixed_pairs[ indexer( Dphi_ij, Deta_ij ) ] += 1.0;
		}
	}

	return;
}
