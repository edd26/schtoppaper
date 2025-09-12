"""
Super, abstract type to bind data sources.
"""
abstract type DataSource end



## ===-===-===-===-===-===-===-===-===-===-
# Matrix loader utils
"""
    PathError

Exception type to handle invalid paths.
"""
struct PathError <: Exception
    msg::String
    invalid_path::String

    function PathError(invalid_path::String)
        new("", invalid_path)
    end

    function PathError(msg::String, invalid_path::String)
        new(msg, invalid_path)
    end
end

Base.showerror(io::IO, e::PathError) = print(io, "Following path could not be found:\n", e.invalid_path, "\n $(e.msg)")


function check_number_ending(name_signature)
    occursin(r"\d+$", name_signature) && "name_signature must not end with a number" |> AssertionError |> throw
end

function check_if_exsits(check_what, check_how::Function)
    MatrixLoaderPathError(path) = PathError("during MatrixLoader creation.", path)
    (check_what |> check_how) || (check_what |> MatrixLoaderPathError |> throw)
end

function get_sample_iterator(data_path, name_signature, file_extension)
    return @pipe readdir(data_path) |>
                 filter(x -> contains(x, name_signature * r"\d+" * ".$(file_extension)" * r"$"), _) |>
                 hcat(split.(_, name_signature)...)[2, :] |>
                 hcat(split.(_, ".$(file_extension)")...)[1, :] |>
                 sort
end
