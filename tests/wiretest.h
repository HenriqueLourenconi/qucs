
#include "schematic_model.h"
#include "place.h"

// dup from lang_qucsator...
static Place const* place_at(pos_t p, Symbol* m)
{
	std::string ps = ":" + std::to_string(getX(p)) + ":" + std::to_string(getY(p));
	auto scope = m->scope();
	assert(scope);
	auto i = scope->find_(ps);
	Place const* ret = nullptr;

	assert(m->mutable_owner());

	if(i == scope->end()){
	}else if(auto p=dynamic_cast<Place const*>(*i)){
		ret = p;
	}else{
		incomplete();// find_again...
		assert(false); //for now.
	}

	if(!ret){
		Symbol* c = symbol_dispatcher.clone("place");
		auto s = prechecked_cast<Place*>(c);
		assert(s);
		s->setPosition(p);
		s->setTypeName("place");
		s->setLabel(ps);
		s->setOwner(m->mutable_owner());
		s->set_port_by_index(0, ps);
		scope->push_back(s);

		ret = s;
	}else{
	}

	return ret;

}

static void connect_push(Element* root, Symbol* sym)
{
	sym->setOwner(root);
// 	M.connect(sym);
	for(unsigned i=0; i<sym->numPorts(); ++i){
		pos_t p = sym->nodePosition(i);
		auto q = place_at(p, sym);
		
		std::string const& l = q->label();
		sym->set_port_by_index(i, l);
	}
	root->scope()->push_back(sym);
}