#ifndef TUPLET_HPP
#define TUPLET_HPP

namespace Bem {

// Tuplet of two numbers a,b which always satisfy a >= b
class Tuplet {
public:
    Tuplet(size_t a_,size_t b_)
        :a(a_),b(b_) { 
            if(a_<b_){
                a = b_;
                b = a_;
            }
        }

    bool operator==(Tuplet const& other) const {
        return (a==other.a and b==other.b); 
    }

    bool operator<(Tuplet const& other) const { // "alphabetical" order
        if(a<other.a) return true;
        if(a==other.a) return b<other.b;
        return false;
    }

    size_t get_a() const { return a; }
    size_t get_b() const { return b; }

private:
    size_t a,b;
};

} // namespace Bem

#endif // TUPLET_HPP