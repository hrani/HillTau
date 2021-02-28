/************************************************************************
 * This program is part of HILLTAU, a framework for fast compact
 * abstractions of signaling events.
 * Copyright	(C) 2021	Upinder S. Bhalla and NCBS
 * It is made available under the terms of the
 * GNU Public License version 3 or later.
 * See the file COPYING.LIB for the full notice.
************************************************************************/

/// #include <pybind11/pybind11.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
using namespace std;
/// namespace py = pybind11;
#include <htHeader.h>
/*
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::map<std::string, double>);
PYBIND11_MAKE_OPAQUE(std::map<std::string, MolInfo>);
PYBIND11_MAKE_OPAQUE(std::map<std::string, ReacInfo>);
PYBIND11_MAKE_OPAQUE(std::map<std::string, EqnInfo>);
*/

MolInfo::MolInfo( const std::string& name_, const std::string& grp_, int order_ = 0, double concInit_ = 0.0 ):
	name(name_),
	grp( grp_ ),
	order( order_),
	concInit( concInit_ ),
	index( 0 )
{
};


ReacInfo::ReacInfo( const string& name_, const string& grp_, 
	const vector< string >& subs_, 
	const map< string, double>& reacObj, 
	const map< string, MolInfo* >& molInfo ):
	name(name_),
	grp( grp_ ),
	KA( reacObj.at("KA") ),
	tau( reacObj.at("tau") ),
	tau2( 1.0 ),
	Kmod( 1.0 ),
	gain( 1.0 ),
	baseline( 0.0 ),
	inhibit( 0 ),
	prdIndex( 0 ),
	subs( subs_ ),
	kh( 1.0 ),
	hillIndex( 0 ),
	reagIndex( 0 ),
	modIndex( 0 ),
	oneSub( false ),
	HillCoeff( 1.0 )
{
	tau2 = tau;
	prdIndex = molInfo.at(name)->index;
	if ( subs.size() == 0 ) {
		throw "Error: Reaction " + name + " has zero reagents\n";
	}
	reagIndex = molInfo.at( subs[0] )->index;
	hillIndex = molInfo.at( subs.back() )->index;
	int numUnique = 1;
	if ( reagIndex != hillIndex ) { // At least two subs
		if (subs.size() == 2) {
			numUnique = 2;
		} else if ( subs.size() > 2 ) {
			if ( subs.back() != subs[1] ) {	// A modifier too
				numUnique = 3;
				modIndex = molInfo.at(subs[1])->index;
			} else { // We have a reagent and multiple Hill ligands, no mod
				numUnique = 2;
			}
		}
	}
	oneSub = ( numUnique == 1);
	HillCoeff = subs.size() + 1 - numUnique;
	kh = pow( KA, HillCoeff);

	auto t = reacObj.find( "tau2" );
	if ( t != reacObj.end() ) {
		tau2 = t->second;
	}
	auto inh = reacObj.find( "inhibit" );
	if ( inh != reacObj.end() ) {
		inhibit = inh->second;
	}
	auto i = reacObj.find( "Kmod" );
	if ( i != reacObj.end() ) {
		if ( numUnique != 3) {
			cout << "Warning: Reaction " << name << " has <3 reagents but Kmod has been specified. Ignoring.\n";
			Kmod = 1.0;
		} else {
			Kmod = i->second;
		}
	} else if ( numUnique == 3)  {
		throw "Error: Reaction " + name + " has 3 reagents but no Kmod specified for modifier.\n";
	}
	auto b = reacObj.find( "baseline" );
	if ( b != reacObj.end() ) {
		baseline = b->second;
	}
	auto g = reacObj.find( "gain" );
	if ( g != reacObj.end() ) {
		gain = g->second;
	}
}

double ReacInfo::eval( Model* model, double dt ) const
{
	double orig = model->conc[ prdIndex ] - baseline;
	double delta = concInf( model->conc ) - orig;
	if ( delta >= 0.0 ) {
		delta *= 1.0 - exp( -dt/tau );
	} else {
		delta *= 1.0 - exp( -dt/tau2 );
	}
	double ret = baseline + orig + delta;
	if (ret < 0.0 ) {
		throw "Error: negative value on: " + name;
	}
	model->conc[ prdIndex ] = ret;
	return ret;
}

double ReacInfo::concInf( const vector< double >& conc ) const
{
	double h = pow( conc[ hillIndex ], HillCoeff );
	double mod = 1.0;
	if ( modIndex != 0 ) {
		mod = conc[ modIndex ] / Kmod;
	}
	if ( oneSub ) {
		return h / KA;
	}

	double s = conc[ reagIndex ] * gain;
	h *= mod;
	if ( inhibit ) {
		return s * (1.0 - h / (h + kh ) );
	} else {
		return s * h / (h + kh);
	}
}

EqnInfo::EqnInfo( const string& name_, const string& grp_, 
			const string& eqnStr_, 
			const map< string, MolInfo* >& molInfo ):
	name(name_),
	grp( grp_ ),
	eqnStr( eqnStr_ ),
	molIndex( 0 )
{
};

double EqnInfo::eval( vector< double >& conc ) const
{
	conc[molIndex] = 0.0;
	return 0.0;
}

string EqnInfo::findMolToken( const string& eqn )
{
	return "";
}

Model::Model()
{;
}

void Model::setReacSeqDepth( int maxDepth )
{
	if ( maxDepth < 1 )
		throw( "Error: maxDepth must be >= 1" );
	sortedReacInfo.clear();
	sortedReacInfo.resize( maxDepth );
}

void Model::assignReacSeq( const string& name, int seq )
{
	auto ri = reacInfo.at( name ); // Assume it is good.
	sortedReacInfo[seq].push_back( ri );
}

void Model::advance( double runtime, int settle )
{
	if (runtime < 1e-6) return;
	double newdt = dt;
	int ratio = 1;
	if (settle) {
		newdt = runtime / 10.0;
		ratio = 10;
	} else {
		if ( dt >= runtime / 2.0 )
			newdt = pow( 10.0, floor( log10( runtime / 2.0 ) ) );
		ratio = int( round( dt / newdt ) );
	}

	int i = 0;
	for (double t = 0.0; t < runtime; t += newdt ) {
		for (auto r = sortedReacInfo.begin(); r != sortedReacInfo.end(); 
						r++) {
			for (auto ri = r->begin(); ri != r->end(); ri++ ) {
				(*ri)->eval( this, newdt );
			}
		}
		for (auto e = eqnInfo.begin(); e != eqnInfo.end(); ++e ) {
			e->second->eval( conc );
		}
		if ( i % ratio == 0 ) {
			plotvec.push_back( conc );
		}
		i++;
	}
	currentTime += runtime;
}

void Model::allocConc()
{
	concInit.resize( molInfo.size(), 0.0 );
	for ( auto m = molInfo.begin(); m != molInfo.end(); m++ ) {
		concInit[ m->second->index ] = m->second->concInit;
	}
	conc = concInit;
}

void Model::reinit()
{
	allocConc();
	for (auto r = sortedReacInfo.begin(); r != sortedReacInfo.end(); r++) {
		for (auto ri = r->begin(); ri != r->end(); ri++) {
			unsigned int j = (*ri)->prdIndex;
			if ((*ri)->inhibit ) {
				concInit[j] = (*ri)->concInf( concInit ) + (*ri)->baseline;
				if ( concInit[j] < 0.0 )
					concInit[j] = 0.0;
			} else {
				concInit[j] = (*ri)->baseline;
			}
		}
	}
	conc = concInit;
	plotvec.clear();
	plotvec.push_back( conc );
}

void Model::makeReac( const string & name, const string & grp, 
				const vector< string >& subs, 
				const map< string, double >& reacObj )
{
	auto r = new ReacInfo( name, grp, subs, reacObj, molInfo );
	reacInfo[ name ] = r;
	if ( molInfo[name]->order == -1 ) {
		if ( r->inhibit ) {
			concInit[ r->prdIndex ] = r->concInf( concInit ) + r->baseline;
			if ( concInit [ r->prdIndex ] < 0.0 )
				concInit [ r->prdIndex ] = 0.0;
		} else {
			concInit[ r->prdIndex ] = r->baseline;
		}
	}
}
void Model::makeMol( const string & name, const string & grp, int order, double concInit )

{
	auto mi = molInfo.find( name );
	if ( mi == molInfo.end() ) { // Make new one.
		auto m = new MolInfo( name, grp, order, concInit );
		m->index = molInfo.size();
		molInfo[ name ] = m;
	} else {
		// Could be the second pass assignment of species with concInits.
		mi->second->concInit = concInit;
		mi->second->order = order;
	}
}

void Model::makeEqn( const string & name, const string & grp, const string& expr )

{
	auto e = new EqnInfo( name, grp, expr, molInfo );
	eqnInfo[ name ] = e;
}

void Model::setConc( unsigned int index, double value )
{
	if ( index > conc.size() )
			throw "Error: index of conc[] out of range";
	conc[index] = value;
}
