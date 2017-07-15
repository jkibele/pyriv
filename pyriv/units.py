import numericalunits as nu
nu.reset_units('SI')

# There may be a more graceful way to do this, but this'll work for now.
# Settings:
proj_unit = 'm'
output_unit = 'km'

in_unit = getattr(nu, proj_unit)
out_unit = getattr(nu, output_unit)

def length_in_display_units(length, in_units=in_unit, out_units=out_unit):
    if in_units.__class__.__name__ == 'str':
        in_units = getattr(nu, in_units)
    if out_units.__class__.__name__ == 'str':
        out_units = getattr(nu, out_units)
    len_out = (length * in_units) / out_units
    return len_out