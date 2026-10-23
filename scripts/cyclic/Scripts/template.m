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
// Function to get code characteristics k1, k2, and dimension
function characteristics(n, generator)
    C := CyclicCode(n, generator);
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
    degreeGen := Degree(generator);
    if n gt 61 then
        if Minimum(n - degreeGen, degreeGen) le 30 then
            MinLeeW := MinimumLeeWeight(C : MaximumTime := TimeLimit);
            if MinLeeW eq -1 then
                skippedTime +:= 1;
            end if;
        else
            skippedDimCap +:= 1;
        end if;
    else
        MinLeeW := MinimumLeeWeight(C); // No time limit when n is smaller or equal to 61
    end if;
	penu_list:=[n,Dim, k1, k2, MinLeeW];
	u_list:=[];
		for item in penu_list do
			Append(~u_list,"\"" cat IntegerToString(item) cat "\"");
		end for;
	gen_string:="\"" cat Sprint(generator) cat "\"" ;
	Append(~u_list,gen_string);
    print u_list;
    return u_list;
end function;

// Divisors function
function GeneratePolynomials(n)
    Z4 := IntegerRing(4);  // Declare Z4 here
    P4<x> := PolynomialRing(Z4);  // Declare P4 here
    Fac := CyclotomicFactors(Z4, n);
    l := #Fac;
    npossibilities := 3^l - 1;
    posslist := [];
    for i in [0..npossibilities] do
        representation := BaseNRepresentation(i, 3, l);
        Append(~posslist, representation);
    end for;
    generator_polynomials := [];
    for j in [1..#posslist] do
        combination := posslist[j];
        function_f := P4!1;
        function_g := P4!1;
        function_h := P4!1;
        gen_poly := P4!1;
        for i in [1..l] do
            if combination[i] eq 0 then
                function_f *:= Fac[i];
            elif combination[i] eq 1 then
                function_g *:= Fac[i];
            elif combination[i] eq 2 then
                function_h *:= Fac[i];
            end if;
        end for;
        gen_poly := (function_f * function_h) + 2*function_f;
        Append(~generator_polynomials, gen_poly);
    end for;
    return generator_polynomials;
end function;

// Procedure to write the header to the file
procedure Opener(file, data_head)
    for header in data_head do
        Write(file, header cat ",");
    end for;
    Write(file, "\n");
end procedure;

// Procedure to write data to the file
procedure Write2File(file, data_rows)
	CSV_String:="";
    for data in data_rows do
	CSV_String:=CSV_String cat data cat ",";
       end for;
	CSV_String := CSV_String[1..#CSV_String-1] cat "\n";
    Write(file, CSV_String);
end procedure;

// Main function
procedure main()
    data_Header := ["n", "Dimension", "k1", "k2", "Distance", "generator"];
    n := {n};
    filename := "/home/belghith1/Exhaustive_Z4/Cyclic/Output/Len" cat IntegerToString(n) cat ".csv";
    
    // Open the file for writing
    file := Open(filename, "w");
    
    Opener(file, data_Header);
    polyList := GeneratePolynomials(n);
    for possgen in polyList do
        code_data := characteristics(n, possgen);
        Write2File(file, code_data);
    end for;
    Flush(file); // Correct way to close the file in Magma
end procedure;

// Run the main procedure
main();
