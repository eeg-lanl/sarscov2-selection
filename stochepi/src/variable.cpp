#include "variable.hpp"


void Variable::make_discrete() {
  if ( !is_discrete ) { // only modify if not already discrete
    disc_value = int(cont_value);
    is_discrete = true;
  }
}

void Variable::make_continuous() {
  if ( is_discrete ) { // only modify if not already continuous
    cont_value = double(disc_value);
    is_discrete = false;
  }
}

double Variable::value() const {
  if ( is_discrete ) {
    // can't use ( ? : ) operator with unequal types int and double
    return disc_value;
  } else {
    return cont_value;
  }
}

void Variable::print(std::ostream & os) const {
  if ( is_discrete ) {
    os << "<var "
       << "type='int' "
       << "val='" << disc_value << "' "
       << "/>";
  } else {
    os << "<var "
       << "type='real' "
       << "val='" << cont_value << "' "
       << "/>";
  }
}

const Variable & Variable::operator=(double x) {
  is_discrete = false;
  cont_value = x;
  disc_value = int(x);
  return *this;
}

const Variable & Variable::operator=(int x) {
  is_discrete = true;
  disc_value = x;
  cont_value = double(x);
  return *this;
}

Variable::operator double() const {
  return value();
}
