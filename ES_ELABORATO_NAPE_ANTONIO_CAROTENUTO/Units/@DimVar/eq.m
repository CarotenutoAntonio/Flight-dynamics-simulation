function result = eq(v1,v2)

if compatible(v1,v2)
    result = v1.value == v2.value;
end

% 2014-05-14/Sartorius: new.