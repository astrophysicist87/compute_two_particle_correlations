#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

using namespace std;

//=================================================================
constexpr size_t Dphi_bins = 100;
constexpr size_t Deta_bins = 100;

constexpr double pi = 3.1415926535897932384626433;
constexpr double twopi = 2.0*pi;
constexpr double Dphimin = -4.0, Dphimax = 4.0;
constexpr double Detamin = -4.0, Detamax = 4.0;
constexpr double Dphi_bw = twopi/(Dphi_bins+1);
constexpr double Deta_bw = (Detamax-Detamin)/(Deta_bins+1);

constexpr size_t n_mix = 10;

vector<vector<double> > signal_pairs, mixed_pairs;


//=================================================================
inline size_t indexer( size_t i, size_t j )
{
	return i * Deta_bins + j;
}

//=================================================================
size_t get_Dphi_bin( double & Dphi )
{
	while ( Dphi < Dphimin ) Dphi += twopi;
	while ( Dphi < Dphimax ) Dphi -= twopi;
	return static_cast<side_t>( (Dphi - Dphimin)/Dphi_bw );
}

//=================================================================
size_t get_Deta_bin( double & Deta )
{
	return static_cast<side_t>( (Deta - Detamin)/Deta_bw );
}

//=================================================================
void read_in_file( string filename, vector<vector<double> > & event );
void get_signal_pairs( const vector<vector<double> > & event );
void get_mixed_pairs( const vector<vector<double> > & event1,
					  const vector<vector<double> > & event2 );

//=================================================================
int main(int argc, char *argv[])
{
	// for choosing random events below
	vector<size_t> event_indices;
	for (size_t iArg = 1; iArg < argc; iArg++) event_indices.push_back( iArg );

	// loop over all files
	for (size_t iArg = 1; iArg < argc; iArg++)
	{
		string filename = argv[iArg];
		vector<vector<double> > event;
		read_in_file( filename, event );

		// update signal pair distribution
		get_signal_pairs( event );

		// choose n_mix other random events to construct background
		std::vector<size_t> mix_events;
		std::sample(event_indices.begin(), event_indices.end(), 
                std::back_inserter(mix_events), n_mix + 1,	// extra event in case
                std::mt19937{std::random_device{}()});		// one is this event

		// loop over randomly chosen events and form background pairs
		size_t mixCount = 0;
		vector<vector<double> > event_to_mix;
		for ( const size_t & mix_event : mix_events )
		{

			if ( mixCount >= n_mix ) break;
			if ( mix_event == iArg ) continue;

			read_in_file( argv[mix_event], event_to_mix );

			get_mixed_pairs( event, event_to_mix );

			mixCount++;

		} // end loop over mixed events

	} // end loop over all events

	// output ratio of signal to background
	// ...do this at home...

	return 0;
}

//==================================================================
void read_in_file(string filename, vector<vector<double> > & event)
{
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
			event.push_back( vector<double>({phi, eta}) );
		}
	}

	infile.close();
	return;
}

//==================================================================
void get_signal_pairs( const vector<vector<double> > & event )
{
	const size_t event_size = event.size();

	for ( size_t i = 0; i < event_size; i++ )
	for ( size_t j = i+1; j < event_size; j++ )
		signal_pairs[indexer( get_Dphi_bin(i[0]-j[0]),
							  get_Deta_bin(i[1]-j[1]) )]++;

	return;
}


//==================================================================
void get_mixed_pairs( const vector<vector<double> > & event1,
					  const vector<vector<double> > & event2 )
{
	const size_t event1_size = event1.size();
	const size_t event2_size = event2.size();

	for ( size_t i = 0; i < event1_size; i++ )
	for ( size_t j = 0; j < event2_size; j++ )
		mixed_pairs[indexer( get_Dphi_bin(i[0]-j[0]),
							 get_Deta_bin(i[1]-j[1]) )]++;
	return;
}