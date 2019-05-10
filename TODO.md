# TODO list

- Make Constraints typed in their bounds !
- Make Polynomial typed in its coefficients values, Exponent in its variables higher type
- make Variable parametric typed (hence essentially a string)

- Work problem manipulation functions and return types: should `get_slacks` return a point ?
- Improve point iteration and access...
- At .sdp export, check that no unhandled space characters are printed (make *all* exponent->string conversions with `format_string`)