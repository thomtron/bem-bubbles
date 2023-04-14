#ifndef TRIPLET_HPP
#define TRIPLET_HPP

namespace Bem {

class Triplet {
public:
    size_t a,b,c;
    Triplet(size_t a,size_t b,size_t c)
        :a(a),b(b),c(c) {}

    Triplet()
        :Triplet(0,0,0) {}

    // optionally: use *(&a+index)
    size_t& operator[](size_t index) {
        switch(index){
            case 0:
                return a;
                break;
            case 1:
                return b;
                break;
            case 2:
                return c;
                break;
            default: 
                throw std::out_of_range ("out of range");
                break;
        }
    }

    const size_t& operator[](size_t index) const {
        switch(index){
            case 0:
                return a;
                break;
            case 1:
                return b;
                break;
            case 2:
                return c;
                break;
            default: 
                throw std::out_of_range ("out of range");
                break;
        }
    }

    void cyclic_reorder(size_t first) {
        if(first == a) return;
        if(first == b) {
            first = a;
            a = b;
            b = c;
            c = first;
            return;
        }
        if(first == c) {
            first = a;
            a = c;
            c = b;
            b = first;
            return;
        }
        throw std::out_of_range("'first' not included!");
    }

    bool operator==(Triplet const& other) const {
        return a == other.a and b == other.b and c == other.c;
    }

    bool operator!=(Triplet const& other) const {
        return not ((*this)==other);
    }
};

} // namespace Bem

#endif // TRIPLET_HPP