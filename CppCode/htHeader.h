/************************************************************************
 * This program is part of HILLTAU, a framework for fast compact
 * abstractions of signaling events.
 * Copyright	(C) 2021	Upinder S. Bhalla and NCBS
 * It is made available under the terms of the
 * GNU Public License version 3 or later.
 * See the file COPYING.LIB for the full notice.
************************************************************************/

class Model;

class MolInfo
{
	public:
			MolInfo( const string& name, const string& grp, int order, double concInit );
			string name;
			string grp;
			int order;
			double concInit;
			unsigned int index;
};	

class ReacInfo
{
	public:
			ReacInfo( const string& name, const string& grp, 
			const vector< string >& subs, 
			const map< string, double>& reacObj, 
			const map< string, MolInfo* >& molInfo );
			string name;
			string grp;
			double KA;
			double tau;
			double tau2;
			double Kmod;
			double gain;
			double baseline;
			int inhibit;
			int prdIndex;
			double concInf( const vector< double >& conc ) const;
			double eval( Model* model, double dt ) const;
			vector< string > subs;

	private:
			double kh;
			unsigned int hillIndex;
			unsigned int reagIndex;
			unsigned int modIndex;
			bool oneSub;
			double HillCoeff;

};

class EqnInfo
{
	public:
			EqnInfo( const string& name, const string& grp, const string& eqnStr, const map< string, MolInfo *>& molInfo );
			string name;
			string grp;
			string eqnStr;
			double eval( vector<double>& conc ) const;
			static string findMolToken( const string& eqn );
	private:
			// Stuff for parser
			vector< unsigned int > subs;
			unsigned int molIndex;

};

class Model
{
	public:
			Model();
			map< string, MolInfo* > molInfo;
			map< string, ReacInfo* > reacInfo;
			map< string, EqnInfo* > eqnInfo;
			double currentTime;
			double dt;
			vector< double > conc;
			vector< double > concInit;
			vector< vector< double > > plotvec;
			
			void makeMol( const string & name, const string & grp, int order, double concInit );
			void makeReac( const string & name, const string & grp, const vector< string >& subs, const map< string, double >& reacObj );
			void makeEqn( const string & name, const string & grp, const string& expr );
			void setReacSeqDepth( int order );
			void assignReacSeq( const string& name, int seq );
			void advance( double runtime, int settle );
			void allocConc();
			void parseEqns();
			void reinit();
	private:
			vector< vector< const ReacInfo* > > sortedReacInfo;
};
