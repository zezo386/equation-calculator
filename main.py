import re
import math
class Term:
    def __init__(self,term):
        self.coefficent = ""
        for i in term:
            if i.isalpha():
                break
            else: 
                self.coefficent += i
        self.variable = term[len(str(self.coefficent)):]
        if self.coefficent in ['','+']:
            self.coefficent = 1
        elif self.coefficent == '-':
            self.coefficent = -1
        else:
            self.coefficent = float(eval(self.coefficent))
        self.variables_dict = self._parse_variables(self.variable)
        self.variable = self._dict_to_variable_string()
        
    def _parse_variables(self, variable_str):

        if not variable_str:
            return {}
        
        variables = {}
        
        pattern = r'([a-zA-Z])(?:\^(\d+))?'
        matches = re.findall(pattern, variable_str)
        
        for var, exp in matches:
            variables[var] = variables.get(var, 0) + (int(exp) if exp else 1)
        
        return variables
    def _dict_to_variable_string(self):
        if not self.variables_dict:
            return ""
        
        parts = []
        for var in sorted(self.variables_dict.keys()):
            exp = self.variables_dict[var]
            if exp == 1:
                parts.append(var)
            else:
                parts.append(f"{var}^{exp}")
        
        return ''.join(parts)

    def __mul__(self, other):
        new_coeff = self.coefficent * other.coefficent
        
        new_variables = self.variables_dict.copy()
        for var, exp in other.variables_dict.items():
            new_variables[var] = new_variables.get(var, 0) + exp
        
        var_string = self._dict_to_variable_string_from_dict(new_variables)
        
        if not var_string and new_coeff != 0:
            return Term(f"{new_coeff}")
        
        return Term(f"{new_coeff}{var_string}")
    def __floordiv__(self, other):
        new_coeff = self.coefficent / other.coefficent
        
        new_variables = self.variables_dict.copy()
        for var, exp in other.variables_dict.items():
            new_variables[var] = new_variables.get(var, 0) - exp
            if new_variables[var] == 0:
                del new_variables[var]
        
        var_string = self._dict_to_variable_string_from_dict(new_variables)
        
        if not var_string and new_coeff != 0:
            return Term(f"{new_coeff}")
        if new_coeff in [1,-1] and var_string:
            return Term(f"{'-' if new_coeff == -1 else ''}{var_string}")
        
        return Term(f"{new_coeff}{var_string}")

    def _dict_to_variable_string_from_dict(self, variables_dict):
        if not variables_dict:
            return ""
        
        parts = []
        for var in sorted(variables_dict.keys()):
            exp = variables_dict[var]
            if exp == 1:
                parts.append(var)
            else:
                parts.append(f"{var}^{exp}")
        
        return ''.join(parts)
    
    def __str__(self):
        return f"{self.coefficent}{self.variable}"
    def __repr__(self):
        return f"{self.coefficent}{self.variable}"
    def __add__(self,other):
        if self.variable != other.variable:
            raise ValueError("Cannot add terms with different variables")
        result = Term(f"{self.coefficent + other.coefficent}{self.variable}")
        return result
    def __sub__(self,other):
        if self.variable != other.variable:
            raise ValueError("Cannot subtract terms with different variables")
        result = Term(f"{self.coefficent - other.coefficent}{self.variable}")
        return result
    



class Polynomial:
    def __init__(self, polynomial):
        cleaned = polynomial.replace(' ', '').replace('-', '+-')
        if cleaned.startswith('+-'):
            cleaned = '-' + cleaned[2:]
        terms_list = [t for t in cleaned.split('+') if t != '']
        self.terms = [Term(t) for t in terms_list]

    def simplify(self):
        d = {}
        for term in self.terms:
            key = frozenset(term.variables_dict.items())
            if key in d:
                d[key] += term.coefficent
            else:
                d[key] = term.coefficent
        
        self.terms = []
        for key, value in d.items():
            if value != 0:
                var_dict = dict(key)
                var_string = ''.join(f"{var}^{exp}" if exp != 1 else var 
                                   for var, exp in sorted(var_dict.items()))
                self.terms.append(Term(f"{value}{var_string}"))
        self.order()
    
    def order(self):
        self.terms.sort(key=lambda x: (-sum(x.variables_dict.values()), x.variable))

    def factorize(self):
        gcf_term = term_gcf(self.terms)

        self.terms.sort(key=lambda x: -sum(x.variables_dict.values()))
        a = self.terms[0].coefficent  # x^2 coefficient
        b = self.terms[1].coefficent if len(self.terms) > 1 else 0  # x coefficient
        c = self.terms[2].coefficent if len(self.terms) > 2 else 0  # constant
        
        
        # Find factors using quadratic formula or factoring method
        factors = self._find_quadratic_factors(a, b, c)
        
        if factors:
            p1, p2 = factors
            product = p1 * p2
            print(f"({p1})({p2})")
            return p1, p2

        new_terms = []
        for term in self.terms:
            new_term = term // gcf_term
            new_terms.append(new_term)
        
        factored_poly = Polynomial("")
        factored_poly.terms = new_terms
        factored_poly.simplify()
        print(f"{gcf_term}( {factored_poly} )")
        return gcf_term, factored_poly

    def _find_quadratic_factors(self, a, b, c):
        
        # Method 1: Factor by grouping (for a=1)
        if a == 1:
            # Look for factors of c that add up to b
            factors_dict = factor(c)
            for factor1, factor2 in factors_dict.items():
                if factor1 + factor2 == b:
                    p1 = Polynomial(f"x + {factor1}")
                    p2 = Polynomial(f"x + {factor2}")
                    return p1, p2
        
        # Method 2: General factoring for a != 1
        else:
            # Look for factors of a*c that add up to b
            product = a * c
            factors_dict = factor(product)
            for factor1, factor2 in factors_dict.items():
                if factor1 + factor2 == b:
                    # Rewrite middle term and factor by grouping
                    poly1 = Polynomial(f"{a}x^2 + {factor1}x")
                    poly2 = Polynomial(f"{factor2}x + {c}")
                    
                    # Factor each group
                    gcf1 = term_gcf(poly1.terms)
                    gcf2 = term_gcf(poly2.terms)
                    
                    factored1_terms = [term // gcf1 for term in poly1.terms]
                    factored2_terms = [term // gcf2 for term in poly2.terms]
                    
                    factored1 = Polynomial("")
                    factored1.terms = factored1_terms
                    factored1.simplify()
                    
                    factored2 = Polynomial("")
                    factored2.terms = factored2_terms
                    factored2.simplify()
                    
                    # Check if the binomial factors are the same
                    if str(factored1) == str(factored2):
                        p1 = Polynomial(f"{gcf1.coefficent}x + {gcf2.coefficent}")
                        p2 = factored1
                        return p1, p2
        
        return None

    def __eq__(self, other):
        """Check if two polynomials are equal"""
        if not isinstance(other, Polynomial):
            return False
        
        self.simplify()
        other.simplify()
        
        if len(self.terms) != len(other.terms):
            return False
        
        # Compare each term
        for term1, term2 in zip(self.terms, other.terms):
            if (abs(term1.coefficent - term2.coefficent) > 1e-10 or 
                term1.variables_dict != term2.variables_dict):
                return False
        
        return True
    def __mul__(self, other):
        result_terms = []
        for term1 in self.terms:
            for term2 in other.terms:
                result_terms.append(term1 * term2)
        
        result_poly = Polynomial("")
        result_poly.terms = result_terms
        result_poly.simplify()
        return result_poly
    
    def __floordiv__(self, other):
        result_terms = []
        for term1 in self.terms:
            for term2 in other.terms:
                result_terms.append(term1 // term2)
        
        result_poly = Polynomial("")
        result_poly.terms = result_terms
        result_poly.simplify()
        return result_poly
    def __truediv__(self, other):
        return self.__floordiv__(other)

    def __str__(self):
        if not self.terms:
            return '0'
        
        result = ''
        for i, term in enumerate(self.terms):
            term_str = str(term)
            if i == 0:
                result = term_str
            else:
                if term.coefficent > 0:
                    result += f" + {term_str}"
                else:
                    result += f" - {term_str[1:] if term_str.startswith('-') else term_str}"
        return result

    def __add__(self, other):
        result = Polynomial(self.__str__() + " " + other.__str__())
        result.simplify()
        return result

    def __sub__(self, other):
        other_neg = Polynomial(other.__str__())
        for term in other_neg.terms:
            term.coefficent *= -1
        return self + other_neg
    

def factor(n):
    factors = {}
    n = int(n)
    if n == 0:
        return {0:0}
    for i in range(-abs(n),abs(n)+1):
        if i == 0:continue
        if n%i == 0:
            factors[i] = n//i
    return factors

def co_gcf(n1,n2):
    n1_factors = factor(n1)
    n2_factors = factor(n2)
    common_factors = set(n1_factors.keys()).intersection(set(n2_factors.keys()))
    gcf = max(common_factors, key=abs) if common_factors else 1
    return gcf

def term_gcf(terms):
    if not terms:
        return Term("1")
    
    gcf_coeff = abs(terms[0].coefficent)
    gcf_vars = terms[0].variables_dict.copy()
    
    for term in terms[1:]:
        gcf_coeff = co_gcf(gcf_coeff, abs(term.coefficent))
        
        for var in list(gcf_vars.keys()):
            if var in term.variables_dict:
                gcf_vars[var] = min(gcf_vars[var], term.variables_dict[var])
            else:
                del gcf_vars[var]
    
    var_string = ''.join(f"{var}^{exp}" if exp != 1 else var 
                         for var, exp in sorted(gcf_vars.items()))
    
    gcf_term = Term(f"{gcf_coeff}{var_string}")
    return gcf_term

class Equation:
    def __init__(self,equation):
        left,right = equation.split('=')
        self.left = Polynomial(left)
        self.right = Polynomial(right)
        self.solve()
    def solve(self):
        if len(self.left.terms) == 1 and len(self.right.terms ) == 1 and list(self.left.terms[0].variables_dict.values())[0] == 1 and self.right.terms[0].variable == '' and self.left.terms[0].coefficent == 1:
            return f"{self.left.terms[0].variable} = {self.right.terms[0].coefficent}"
        if len(self.right.terms) == 1 and len(self.left.terms ) == 1 and list(self.left.terms[0].variables_dict.values())[0] == 1:
            right_term = self.right.terms[0]
            left_term = self.left.terms[0]
            if right_term.variable == '':
                print("steps:")
                print(f"{variable_term.variable} = {right_term.coefficent}/{left_term.coefficent}")
                print(f"{variable_term.variable} = {right_term.coefficent/left_term.coefficent}")
                self.solution = f"{variable_term.variable} = {right_term.coefficent/left_term.coefficent}"
            return Equation(f"{self.solution}")
        if len(self.right.terms) == 1 and len(self.left.terms ) == 1 and not list(self.left.terms[0].variables_dict.values())[0] == 1:
            exp = list(self.left.terms[0].variables_dict.values())[0]
            right_term = self.right.terms[0]
            left_term = self.left.terms[0]
            variable_term = list(self.left.terms[0].variables_dict.keys())[0]
            if right_term.variable == '':
                print("steps:")
                print(f"{variable_term}^{exp} = {right_term.coefficent}/{left_term.coefficent}")
                print(f"{variable_term} = ({right_term.coefficent/left_term.coefficent})^(1/{exp})")
                if exp%2 == 0 and (right_term.coefficent/left_term.coefficent)<0:
                    print("No real solution exists since we cannot take an even root of a negative number.")
                    self.solution = "No real solution"
                    return self.solution
                if exp%2 == 0:
                    print(f"{variable_term} = ±{round((right_term.coefficent/left_term.coefficent)**(1/exp),5)}")
                    self.solution = f"{variable_term} = {round((right_term.coefficent/left_term.coefficent)**(1/exp),5)}"
                    self.solution2 = f"{variable_term} = -{round((right_term.coefficent/left_term.coefficent)**(1/exp),5)}"
                    print(self.solution)
                    print(self.solution2)
                    return set({Equation(f"{self.solution}"),Equation(f"{self.solution2}")})
                print(f"{variable_term} = {round((right_term.coefficent/left_term.coefficent)**(1/exp),5)}")
                self.solution = f"{variable_term} = {round((right_term.coefficent/left_term.coefficent)**(1/exp),5)}"
                print(self.solution)
                Equation(f"{self.solution}")
        elif len(self.right.terms) == 1 and len(self.left.terms) == 2:
            if self.left.terms[1].variable != '':
                self.left.terms.append(Term(f"{-1*self.right.terms[0].coefficent}"))
                self.right.terms[0].coefficent = 0
                Equation(f"{self.solution}")
            if self.right.terms[0].variable != '':
                self.right.terms.append(Term(f"{-1*self.left.terms[0].coefficent}"))
                self.left.terms.append(Term(f"{-1*self.left.terms[0].coefficent}"))
                
                return Equation(f"{self.solution}")
            right_terms = self.right.terms[0]
            left_terms = self.left.terms
            constant_term = [t for t in left_terms if t.variable == ''][0]
            variable_term = [t for t in left_terms if t.variable != ''][0]
            if right_terms.variable == '':
                print("steps:")
                print(f"{variable_term.coefficent}{variable_term.variable} = {right_terms.coefficent} - {constant_term.coefficent}  ")
                print(f"{variable_term.coefficent}{variable_term.variable} = {right_terms.coefficent-constant_term.coefficent}  ")
                print(f"{variable_term.variable} = {(right_terms.coefficent - constant_term.coefficent)}/{variable_term.coefficent}  ")
                print(f"{variable_term.variable} = {(right_terms.coefficent - constant_term.coefficent)/variable_term.coefficent}  ")
                self.solution = f"{variable_term.variable} = {(right_terms.coefficent - constant_term.coefficent)/variable_term.coefficent}"
                return Equation(f"{self.solution}")
        elif len(self.right.terms) == 2 and len(self.left.terms) == 2:
            right_terms = self.right.terms
            left_terms = self.left.terms
            constant_term_left = [t for t in left_terms if t.variable == ''][0]
            constant_term_right = [t for t in right_terms if t.variable == ''][0]
            variable_term_left = [t for t in left_terms if t.variable != ''][0]
            variable_term_right = [t for t in right_terms if t.variable != ''][0]
            if variable_term_left.variable == variable_term_right.variable:
                print("steps:")
                print(f"{variable_term_left.coefficent}{variable_term.variable} - {variable_term_right.coefficent}{variable_term.variable} = {constant_term_right.coefficent} - {constant_term_left.coefficent}  ")
                print(f"{variable_term_left.coefficent - variable_term_right.coefficent}{variable_term.variable} = {constant_term_right.coefficent - constant_term_left.coefficent}  ")
                print(f"{variable_term.variable} = {(constant_term_right.coefficent - constant_term_left.coefficent)}/{(variable_term_left.coefficent - variable_term_right.coefficent)}  ")
                print(f"{variable_term.variable} = {(constant_term_right.coefficent - constant_term_left.coefficent)/(variable_term_left.coefficent - variable_term_right.coefficent)}  ")
                self.solution = f"{variable_term.variable} = {(constant_term_right.coefficent - constant_term_left.coefficent)/(variable_term_left.coefficent - variable_term_right.coefficent)}"
                return Equation(f"{self.solution}")
        elif len(self.left.terms) == 3:
            if len(self.right.terms) == 1:
                if self.right.terms[0].coefficent == 0 and self.right.terms[0].variable == '':
                    left_poly = Polynomial("")
                    left_poly.terms = self.left.terms
                    left_poly.simplify()
                    left_poly_str = str(left_poly)
                    print("steps:")
                    print(f"{left_poly_str} = 0")
                    print(f"{left_poly_str} factored:")
                    factors = left_poly.factorize()
                    print("Set each factor to zero and solve for x.")
                    self.solution_set = set({Equation(f"{str(factor)} = 0").solution for factor in factors if isinstance(factor, Polynomial)})
                    print(f"final solution set: {self.solution_set}")
                    return self.solution_set
                elif self.right.terms[0].coefficent != 0 and self.right.terms[0].variable == '':
                    left_poly = Polynomial("")
                    left_poly.terms = self.left.terms
                    left_poly.simplify()
                    right_term = self.right.terms[0]
                    left_poly_str = str(left_poly)
                    print("steps:")
                    print(f"{left_poly_str} = {right_term.coefficent}")
                    print(f"{left_poly_str} - {right_term.coefficent} = 0")
                    new_left_poly = Polynomial(f"{left_poly_str} - {right_term.coefficent}")
                    new_left_poly.simplify()
                    new_left_poly_str = str(new_left_poly)
                    print(f"{new_left_poly_str} factored:")
                    factors = new_left_poly.factorize()
                    print("Set each factor to zero and solve for x.")
                    self.solution_set = set({Equation(f"{str(factor)} = 0").solution for factor in factors if isinstance(factor, Polynomial)})
                    print(f"final solution set: {self.solution_set}")
                    return self.solution_set
                elif self.right.terms[0].variable != '':
                    self.right.terms[0].coefficent *= -1
                    left_poly = Polynomial("")
                    left_poly.terms = self.left.terms + [self.right.terms[0]]
                    left_poly.simplify()
                    left_poly_str = str(left_poly)
                    print("steps:")
                    print(f"{left_poly_str} = 0")
                    print(f"{left_poly_str} factored:")
                    factors = left_poly.factorize()
                    print("Set each factor to zero and solve for x.")
                    self.solution_set = set({Equation(f"{str(factor)} = 0").solution for factor in factors if isinstance(factor, Polynomial)})
                    print(f"final solution set: {self.solution_set}")
                    return self.solution_set
            elif len(self.right.terms) == 2:
                self.right.terms[0].coefficent *= -1
                self.right.terms[1].coefficent *= -1
                left_poly = Polynomial("")
                left_poly.terms = self.left.terms + self.right.terms
                left_poly.simplify()
                left_poly_str = str(left_poly)
                print("steps:")
                print(f"{left_poly_str} = 0")
                print(f"{left_poly_str} factored:")
                factors = left_poly.factorize()
                print("Set each factor to zero and solve for x.")
                self.solution_set = set({Equation(f"{str(factor)} = 0").solution for factor in factors if isinstance(factor, Polynomial)})
                print(f"final solution set: {self.solution_set}")
                return self.solution_set
            elif len(self.right.terms) == 3:
                self.right.terms[0].coefficent *= -1
                self.right.terms[1].coefficent *= -1
                self.right.terms[2].coefficent *= -1
                left_poly = Polynomial("")
                left_poly.terms = self.left.terms + self.right.terms
                left_poly.simplify()
                left_poly_str = str(left_poly)
                print("steps:")
                print(f"{left_poly_str} = 0")
                print(f"{left_poly_str} factored:")
                factors = left_poly.factorize()
                print("Set each factor to zero and solve for x.")
                self.solution_set = set({Equation(f"{str(factor)} = 0").solution for factor in factors if isinstance(factor, Polynomial)})
                print(f"final solution set: {self.solution_set}")
                return self.solution_set 
    def __str__(self):
        return f"{self.left} = {self.right}"


e1 = Equation("3x^2 -1 = -13")


