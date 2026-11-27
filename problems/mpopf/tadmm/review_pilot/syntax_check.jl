function reject_parse_errors(ex)
    ex isa Expr || return
    ex.head in (:error,:incomplete) && error(ex)
    foreach(reject_parse_errors,ex.args)
end
for file in ARGS
    reject_parse_errors(Meta.parseall(read(file,String);filename=file))
    println("SYNTAX_OK ",file)
end
