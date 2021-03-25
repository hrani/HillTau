#include <string>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <exprtk.hpp>
using namespace std;

namespace py = pybind11;
#include <htHeader.h>

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::map<string, MolInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, ReacInfo>);
PYBIND11_MAKE_OPAQUE(std::map<string, EqnInfo>);

PYBIND11_MODULE(ht, m) {
	py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::class_<MolInfo>(m, "MolInfo")
        .def( 
			py::init<const std::string &, const std::string &, int, double>(), py::arg("name"), py::arg("grp"), py::arg("order") = 0, py::arg("concInit") = 0.0)
		.def_readwrite("name", &MolInfo::name)
		.def_readwrite("grp", &MolInfo::grp)
		.def_readwrite("order", &MolInfo::order)
		.def_readwrite("concInit", &MolInfo::concInit)
		.def_readwrite("index", &MolInfo::index);

	/////////////////////////////////////////////////////////////////////
    py::class_<ReacInfo>(m, "ReacInfo")
        .def( 
			py::init<const std::string &, const std::string &, const vector< string >&, const map< string, double>&, const map< string, MolInfo*>&>())
		.def_readwrite("name", &ReacInfo::name)
		.def_readwrite("grp", &ReacInfo::grp)
		.def_readwrite("KA", &ReacInfo::KA)
		.def_readwrite("tau", &ReacInfo::tau)
		.def_readwrite("tau2", &ReacInfo::tau2)
		.def_readwrite("Kmod", &ReacInfo::Kmod)
		.def_readwrite("gain", &ReacInfo::gain)
		.def_readwrite("baseline", &ReacInfo::baseline)
		.def_readwrite("inhibit", &ReacInfo::inhibit)
		.def_readwrite("prdIndex", &ReacInfo::prdIndex)
		.def_readwrite("kh", &ReacInfo::kh)
		.def_readonly("HillCoeff", &ReacInfo::HillCoeff)
		.def_readonly("subs", &ReacInfo::subs)
		.def( "eval", &ReacInfo::eval, "Evaluator for Reacs" )
		.def( "concInf", &ReacInfo::concInf, "Computes steady-state value of reaction output" );

	/////////////////////////////////////////////////////////////////////
    py::class_<EqnInfo>(m, "EqnInfo")
        .def( 
			py::init<const std::string &, const std::string &, const std::string&, const vector< string >&, const map< string, MolInfo* >&, vector< double >&>())
		.def_readwrite("name", &EqnInfo::name)
		.def_readwrite("grp", &EqnInfo::grp)
		.def_readwrite("eqnStr", &EqnInfo::eqnStr)
		.def_readonly("subs", &EqnInfo::subs)
		.def( "eval", &EqnInfo::eval, "Evaluator for Eqns" );
	/////////////////////////////////////////////////////////////////////

    py::class_<Model>(m, "Model")
        .def(py::init())
		.def_readwrite("molInfo", &Model::molInfo)
		.def_readwrite("reacInfo", &Model::reacInfo)
		.def_readwrite("eqnInfo", &Model::eqnInfo)
		.def_readonly("currentTime", &Model::currentTime)
		.def_readwrite("dt", &Model::dt)
		.def_readwrite("conc", &Model::conc)
		.def_readwrite("concInit", &Model::concInit)
		.def_readonly("plotvec", &Model::plotvec)
		.def( "makeMol", &Model::makeMol, "Create MolInfo object.", py::arg("name"), py::arg("grp"), py::arg("order") = 0, py::arg("concInit") = 0.0, py::arg("preserveConcInit") = false )
		.def( "makeReac", &Model::makeReac, "Create ReacInfo object.", py::arg("name"), py::arg("grp"), py::arg("subs"), py::arg("reacParms"))
		.def( "makeEqn", &Model::makeEqn, "Create EqnInfo object.", py::arg("name"), py::arg("grp"), py::arg("expr"), py::arg( "eqnSubs" ) )
		.def( "setReacSeqDepth", &Model::setReacSeqDepth, "Defines how deep is the sequence of reactions, that is, the size of sortedReacInfo.")
		.def( "assignReacSeq", &Model::assignReacSeq, "Builds up sortedReacOrder vectors.")
		.def( "advance", &Model::advance, "Advances the simulation", py::arg( "runtime" ), py::arg( "settle" ) = 0 )
		.def( "reinit", &Model::reinit, "Reinits all conc values" )
		.def( "allocConc", &Model::allocConc, "Allocates and initializes conc vectors" );
}

