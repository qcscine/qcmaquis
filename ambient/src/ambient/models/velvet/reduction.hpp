namespace ambient { namespace models { namespace velvet {

    inline reduction::reductionq::reductionq(){

    }

    inline void reduction::reductionq::push(layout::entry* e){

    }

    inline reduction::reduction(revision* r)
    : target(r)
    {
        dim2 dim = this->target->get_layout().get_grid_dim();
        for(size_t x = 0; x < dim.x; x++){
            entries.push_back(std::vector<reductionq*>());
            for(size_t y = 0; y < dim.y; y++){
                entries[x].push_back(NULL);
            }
        }
    };

    inline reduction::~reduction(){
    }

    inline layout::entry& reduction::operator()(size_t x, size_t y){
        if(this->entries[x][y] == NULL) this->entries[x][y] = new reductionq();
        layout::entry* block = NULL; //controller.alloc_block(*this->target);
        this->entries[x][y]->push(block);
        return *block;
    }

} } }

