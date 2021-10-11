#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

//=================================================================
constexpr size_t Dphi_bins = 32;
constexpr size_t Deta_bins = 28;
constexpr size_t n_mix = 100;

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
void read_in_file( string filename, vector<double> & event );
void get_signal_pairs( const vector<double> & event );
void get_mixed_pairs( const vector<double> & event1, const vector<double> & event2 );

//=================================================================
int main(int argc, char *argv[])
{
	signal_pairs.resize(Dphi_bins*Deta_bins);
	mixed_pairs.resize(Dphi_bins*Deta_bins);

	// for choosing random events below
	vector<size_t> event_indices;
	for (size_t iArg = 1; iArg < argc; iArg++) event_indices.push_back( iArg );

	// loop over all files
	for (size_t iArg = 1; iArg < argc; iArg++)
	{
		cout << "Reading in " << argv[iArg] << "; " << flush;

		string filename = argv[iArg];
		vector<double> event;
		read_in_file( filename, event );

		cout << "read in " << event.size()/2 << " particles." << endl;

		// update signal pair distribution
		cout << "\t - updating signal distribution..." << flush;
		get_signal_pairs( event );
		cout << "done!" << endl;

		// choose n_mix other random events to construct background
		std::vector<size_t> mix_events;
		std::sample(event_indices.begin(), event_indices.end(), 
                std::back_inserter(mix_events), n_mix + 1,	// extra event in case
                std::mt19937{std::random_device{}()});		// one is this event

		// loop over randomly chosen events and form background pairs
		size_t mixCount = 0;
		vector<double> event_to_mix;
		for ( const size_t & mix_event : mix_events )
		{
			if ( mixCount >= n_mix ) break;
			if ( mix_event == iArg ) continue;

			cout << "\t - mixing " << argv[iArg] << " with " << argv[mix_event] << endl;

			read_in_file( argv[mix_event], event_to_mix );

			cout << "\t - updating mixed distribution..." << flush;
			get_mixed_pairs( event, event_to_mix );
			cout << "done!" << endl;

			mixCount++;

		} // end loop over mixed events

	} // end loop over all events

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
void read_in_file(string filename, vector<double> & event)
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

			// impose cuts on eta
			if ( abs(eta) > 2.4 ) continue;
			event.push_back(phi);
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
		const size_t Dphi_ij = get_Dphi_bin(event[2*i+0]-event[2*j+0]);
		const size_t Deta_ij = get_Deta_bin(event[2*i+1]-event[2*j+1]);
		if ( Dphi_ij >= 0 && Dphi_ij < Dphi_bins && Deta_ij >= 0 && Deta_ij < Deta_bins )
			signal_pairs[ indexer( Dphi_ij, Deta_ij ) ] += 1.0;

		// ensure reflection symmetry under eta --> -eta (is this correct?)
		const size_t Dphi_ji = get_Dphi_bin(event[2*j+0]-event[2*i+0]);
		const size_t Deta_ji = get_Deta_bin(event[2*j+1]-event[2*i+1]);
		if ( Dphi_ji >= 0 && Dphi_ji < Dphi_bins && Deta_ji >= 0 && Deta_ji < Deta_bins )
			signal_pairs[ indexer( Dphi_ji, Deta_ji ) ] += 1.0;
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
		const size_t Dphi_ij = get_Dphi_bin(event1[2*i+0]-event2[2*j+0]);
		const size_t Deta_ij = get_Deta_bin(event1[2*i+1]-event2[2*j+1]);
		if ( Deta_ij >= 0 && Deta_ij < Deta_bins )
			mixed_pairs[ indexer( Dphi_ij, Deta_ij ) ] += 1.0;
	}

	return;
}
