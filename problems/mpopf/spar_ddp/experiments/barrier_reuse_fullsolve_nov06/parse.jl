using Test
function errors!(ex,path)
 if ex isa Expr
  @assert !(ex.head in (:error,:incomplete)) "Parse error $path $ex"
  for x in ex.args;errors!(x,path);end
 end
end
for p in ARGS
 if isdir(p)
  for (d,_,fs) in walkdir(p),f in fs
   endswith(f,".jl")||continue;file=joinpath(d,f);errors!(Meta.parseall(read(file,String)),file)
  end
 else
  errors!(Meta.parseall(read(p,String)),p)
 end
end
println("BR_PARSE_PASS")
