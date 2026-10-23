///////////
function BaseNRepresentation(number, base, length)
    digits := [];
    quotient := number;

    for i in [1..length] do
        // Perform the division and store the remainder
        remainder := quotient mod base;
        quotient := quotient div base;
        Append(~digits, remainder);
    end for;

    // Reverse the digits to get the correct order
    digits := Reverse(digits);
    
    return digits;
end function;


function polynomials_to_coeffs(factors_List)
    Coeff_List := [];
    for poly in factors_List do
        Append(~Coeff_List, Coefficients(poly));
    end for;
    return Coeff_List;
end function;

function Monic_b_bar_changes(coeff_list)
    changes_list := [];
    num_combinations := 2^(#coeff_list);
    for i in [0..num_combinations-1] do
        Append(~changes_list, BaseNRepresentation(i, 2, #coeff_list));
    end for;
    // Filter out elements that start with 1
    filtered_list := [x : x in changes_list | x[#coeff_list] eq 0];
    return filtered_list;
end function;

function b_bar_lister(coeff_list)
    Z4 := IntegerRing(4);  // Declare Z4 here
    P4<x> := PolynomialRing(Z4);  // Declare P4 here
    // we take b_bar to be the first one in Coeff list.
    possibilities_list := Monic_b_bar_changes(coeff_list);
    possible_b_bars := [];
    for poss in possibilities_list do
        b_bar_nm := coeff_list;
        for i in [1..#b_bar_nm] do
            if poss[i] eq 1 then
                b_bar_nm[i] := Z4!(Integers()!b_bar_nm[i] + 2);  // Add 2 and convert back to Z4
            end if;
        end for;
        Append(~possible_b_bars, b_bar_nm);
    end for;
    return possible_b_bars;
end function;

function Co2poly(coeff_list, power_indice);
    Z4 := IntegerRing(4);  // Declare Z4 here
    P4<x> := PolynomialRing(Z4);  // Declare P4 here
    l := #coeff_list;
    polynomial_needed := P4!0;
    for i in [1..l] do
        poly_to_add := coeff_list[i] * x^((i-1) * power_indice);
        polynomial_needed := polynomial_needed + poly_to_add;
    end for;
    return polynomial_needed;
end function;



// Function to pad the coefficients of a polynomial to a specified length
function PadCoefficients(coeffs, length);
    padded_coeffs := coeffs;
    while #padded_coeffs lt length do
        Append(~padded_coeffs, 0);
    end while;
    return padded_coeffs;
end function;

// Function to create a circulant matrix from a polynomial
function Circ(poly, num_columns);
Z4 := IntegerRing(4);  // Declare Z4 here
    P4<x> := PolynomialRing(Z4);  // Declare P4 here
    coeffs := Coefficients(poly);
    coeffs := PadCoefficients(coeffs, num_columns);
    n := #coeffs;
    C := ZeroMatrix(Parent(coeffs[1]), n, n);
    for i in [1..n] do
        for j in [1..n] do
            C[i][j] := coeffs[(j - i) mod n + 1];
        end for;
    end for;
    return C;
end function;

function give_gen_Matrix(generator_list,code_length);
    Z4 := IntegerRing(4);  // Declare Z4 here
    P4<x> := PolynomialRing(Z4);  // Declare P4 here
    G_1:=Circ(generator_list[1],code_length);
    G_2:=Circ(generator_list[2],code_length);
    big_G:=VerticalJoin(G_1,G_2);
return big_G;
end function;
    
function find_and_remove_sublist(ultimate_list, k1, k2)
    index := 1;
    while index le #ultimate_list do
        sublist := ultimate_list[index];
        if sublist[2] eq k1 and sublist[3] eq k2 then
            removed_list := sublist;
            Remove(~ultimate_list, index);
            return removed_list, ultimate_list;
        end if;
        index +:= 1;
    end while;
    return [], ultimate_list;  // If no matching sublist is found
end function;

function add_code_to_ultimate_list(ultimate_list, n, k1, k2, weight, polynomials)
    sublist, updated_ultimate_list := find_and_remove_sublist(ultimate_list, k1, k2);
    
    if #sublist eq 0 then
        // If no matching sublist is found, create a new one
        new_entry := <n, k1, k2, <weight, [polynomials]>,[]>;
        Append(~updated_ultimate_list, new_entry);
    else
        existing_weight:=sublist[4][1];
        if weight eq -1 then
                list_minus1_codes_generators:=sublist[5];
                Append(~list_minus1_codes_generators, [polynomials]);
                sublist:=<sublist[1],sublist[2].sublist[3],sublist[4],list_minus1_codes_generators>;
        elif weight gt existing_weight then
            sublist:=<sublist[1],sublist[2],sublist[3], <weight, [polynomials]>,sublist[5]>;
        elif weight eq existing_weight then
            Append(~sublist[4][2],polynomials);
        end if;
        // Re-append the modified sublist to the updated_ultimate_list
        Append(~updated_ultimate_list, sublist);
    end if;

    return updated_ultimate_list;
end function;

function characteristics(C, polynomials, ultimate_list)
    TimeLimit := 10;
    skippedTime := 0;
    skippedDimCap := 0;
    len := Length(C);
    k1 := 0;
    B := Basis(C);
    for i in [1..#B] do
        j := 1;
        while B[i][j] eq 0 do
            j := j + 1;
        end while;
        if B[i][j] eq 1 then
            k1 := k1 + 1;
        end if;
    end for;
    Dim := PseudoDimension(C);
    k2 := Dim - k1;
    if len gt 61 then
        if Minimum(len - Dim, Dim) le 30 then
            MinLeeW := MinimumLeeWeight(C : MaximumTime := TimeLimit);
            if MinLeeW eq -1 then
                skippedTime +:= 1;
            end if;
        else
            skippedDimCap +:= 1;
        end if;
    else
        MinLeeW := MinimumLeeWeight(C); // No time limit when len is smaller or equal to 61
    end if;

    ultimate_list := add_code_to_ultimate_list(ultimate_list, len, k1, k2, MinLeeW, polynomials);

    return ultimate_list;
end function;


function Finisher(n,genlist,ultimate_list);
    Big_G:=give_gen_Matrix(genlist, 2*n);
    new_ultimate_list:=characteristics(LinearCode(Big_G,genlist,ultimate_list));
return new_ultimate_list;



//code to fix the problem
//for item in underultimate_list do
			//Append(~u_list,"\"" cat IntegerToString(item) cat "\"");
		//end for;
	//gen_string:="\"" cat Sprint(generator) cat "\"" ;
//u_list:=[];
		
	//Append(~u_list,gen_string);
    //print u_list;

function Computation(n)
    Z4 := IntegerRing(4);
    P4<x> := PolynomialRing(Z4);
    Original_Factors := CyclotomicFactors(Z4, n);
    r := #Original_Factors;
    ultimate_list:= [];
    
    for i_0 in [0..(2^r)-1] do
        d_factors_absorption := BaseNRepresentation(i_0, 2, r);
        d_factors_list := [];
        for j in [1..r] do
            if d_factors_absorption[j] eq 1 then 
                Append(~d_factors_list, Original_Factors[j]);
            end if;
        end for;
        
        Dless_Factors := [x : x in Original_Factors | x notin d_factors_list];
        r_1 := #Dless_Factors;
        print "d has absorped the factors:",d_factors_list;

        print "New D-Level initiated.";

        if r_1 eq 0 then
            continue;
        end if;
        
        
        for i_1 in [0..(2^(r_1)-1)] do
        print "New B-Level Initiated.";
        b:=P4!1;
            b_factors_absorption := [];
            b_factors_absorption := BaseNRepresentation(i_1, 2, r_1);
            b_factors_list := [];
            for j in [1..r_1] do
                if b_factors_absorption[j] eq 1 then 
                    Append(~b_factors_list, Dless_Factors[j]);
                    b *:= Dless_Factors[j];
                end if;
            end for;
            Remaining_Factors := [x : x in Dless_Factors | x notin b_factors_list];
            r_2:=#Remaining_Factors;
            All_it_b_tilde:=b_bar_lister(Coefficients(b));
            // if no factors are left
            if r_2 eq 0 then
                for b_tilde in All_it_b_tilde do
                    trivial_2:=[];
                    Append(~trivial_2, Co2poly(b_tilde, 1));
                    Append(~trivial_2, 2*Co2poly(Coefficients(b),1));
                    ultimate_list:=Finisher(n,trivial_2,ultimate_list);
                end for;
            continue;
            end if;
            final_factors_representation:=[];
            for i_3 in [0..4^(r_2) - 1] do
                single_list:=BaseNRepresentation(i_3,4,r_2);
                factor_list:=[];
                a_1:=P4!1;
                a_2:=P4!1;
                a_3:=P4!1;
                c:=P4!1;
                for j in [1..#single_list] do
                    if single_list[j] eq 0 then
                        a_1*:=Remaining_Factors[j];
                        elif single_list[j] eq 1 then
                        a_2*:=Remaining_Factors[j];
                        elif single_list[j] eq 2 then
                        a_3*:=Remaining_Factors[j];
                        elif single_list[j] eq 3 then
                        c*:=Remaining_Factors[j];
                    end if;
                end for;
                Append(~factor_list,Coefficients(a_1));
                Append(~factor_list,Coefficients(a_2));
                Append(~factor_list,Coefficients(a_3));
                Append(~factor_list,Coefficients(c));
                Append(~final_factors_representation,factor_list);
            end for;
            for b_tilde_coeff in All_it_b_tilde do
                for factor_distribution in final_factors_representation do
                genlist:=[];
                print "The components for this iteration are: ";
                print "d absorbed", d_factors_list;
                print "b absorbed", b_factors_list;
                print "b is ", b;
                print "the b tilde of this iteration is: ", Co2poly(b_tilde_coeff,1);
                print "the remaining factors are ", Remaining_Factors;
                print "a_1 is ", Co2poly(factor_distribution[1],1);
                print "a_2 is ", Co2poly(factor_distribution[2],1);
                print "a_3 is ", Co2poly(factor_distribution[3],1);
                print "c is ", Co2poly(factor_distribution[4],1);
                f_1:=P4!1;
                f_2:=P4!1;
                f_1*:=Co2poly(factor_distribution[1],2)*Co2poly(factor_distribution[2],2)*Co2poly(factor_distribution[3],2)*Co2poly(b_tilde_coeff,1)*Co2poly(factor_distribution[4],1) mod (x^(2*n) - 1);
                f_2*:=2*Co2poly(factor_distribution[1],2)*Co2poly(factor_distribution[2],1)*b mod (x^(2*n) - 1);
                Append(~genlist,f_1);
                Append(~genlist,f_2);
                ultimate_list:=Finisher(n,genlist,ultimate_list);
                end for;
            end for;
        end for;
    
    end for;
    return ultimate_list;
end function;

// Procedure to write the header to the file
procedure Opener(file, data_head)
    for header in data_head do
        Write(file, header cat ",");
    end for;
    Write(file, "\n");
end procedure;

procedure Write2File(file, n, k1, k2, distance, generators)
    CSV_String := IntegerToString(n) cat "," cat IntegerToString(k1) cat "," cat IntegerToString(k2) cat "," cat distance cat ",";
    generator_strings := [];
    for generator in generators do
        generator_str := "[" cat Join([Sprint(gen) : gen in generator], ", ") cat "]";
        Append(~generator_strings, generator_str);
    end for;
    CSV_String := CSV_String cat "[" cat Join(generator_strings, ", ") cat "]\n";
    Write(file, CSV_String);
end procedure;



procedure main()
    data_Header := ["n", "k1", "k2", "Distance", "generator"];
    n := 5;  // Example value for n, replace as needed
    filename := "/home/belghith1/Exhaustive_Z4/Cyclic/OddlyEven/Output/Length" cat IntegerToString(n) cat ".csv";
    
    // Open the file for writing
    file := Open(filename, "w");
    
    // Write the header
    Opener(file, data_Header);
    
    code_list := Computation(n);
    for code in code_list do
        n := code[1];
        k1 := code[2];
        k2 := code[3];
        distance := IntegerToString(code[4][1]);
        generators := code[4][2];
        
        // Write the main part of the codes
        for generator in generators do
            Write2File(file, n, k1, k2, distance, [generator]);
        end for;
        
        // Write the -1 weight codes
        distance := "-";
        minus1_generators := code[5];
        for generator in minus1_generators do
            Write2File(file, n, k1, k2, distance, generator);
        end for;
    end for;
    
    Flush(file);
end procedure;

main();
