'''
Define a bunch a accepted units for several physical quantities (temperature, ...)
'''


class UnknownUnitError(Exception):
    "Raised when units is not reckognized among a set of units"
    pass



class units():
    '''
    "units" object contains a reference unit name, and a list of accepted units name,
    and corresponding list factors and offsets to convert the accepted units in the
    reference unit.
    '''

    def __init__(self, reference_units, accepted_units, conversion_factor, conversion_offset=None):
        '''
        Create a "units" object.
        Mandatory arguments:
          * reference_units: string. Name of the reference units
          * accepted_units: list. Names of the accepted units. The element of the list
                            must be strings, of lists of string (if several accepted
                            units have the same conversion parameters).
          * conversion_factor: list of numbers (real of integer): factor to convert the accepted
                               unit into the reference unit
        Optional arguments:
          * conversion_offset: list of numbers (real of integer): offset to convert the accepted
                               unit into the reference unit. By default, it is 0.
        
        NOTE: the unit conversion is 'factor*value + offset'
        '''

        # Input arguments
        # ---------------


        # reference units

        if isinstance(reference_units, str):
            self.reference = reference_units
        else:
            raise TypeError('"reference_units" argument must be a string')


        # list of accepted units

        if isinstance(accepted_units, list):
            passed = True
            for element in accepted_units:
                if isinstance(element, str):
                    loc_passed = True
                elif isinstance(element, list):
                    loc_passed = True
                    for el in element:
                        if not isinstance(el, str):
                            loc_passed = False

                else:
                    loc_passed = False

                if not loc_passed:
                    passed = False

        else:
            passed = False
                
        if passed:
            self.accepted = accepted_units
        else:
            raise TypeError('"accepted_units" argument must be a list of strings, or of lists of strings')


        # List of conversion factors

        if isinstance(conversion_factor, list):
            if len(conversion_factor) == len(accepted_units):
                self.factor = conversion_factor
            else:
                raise ValueError('"accepted_units" and "conversion_factor" arguments must have the same length')
        
        else:
            raise TypeError('"conversion_factor" must be a list of numbers')


        # List of conversion offsets

        if conversion_offset is None:

            self.offset = len(accepted_units)*[0]

        else:

            if isinstance(conversion_offset, list):
                if len(conversion_offset) == len(accepted_units):
                    self.offset = conversion_offset
                else:
                    raise ValueError('"accepted_units" and "conversion_offset" arguments must have the same length')
            
            else:
                raise TypeError('"conversion_offset" must be a list of numbers')




    # Class methods
    # -------------

    def convert(self, value, unit):
        '''
        convert the value "value" of a given unit "unit" into the reference unit
        of the object, if it belongs tothe list of accepted units.
        "value" must be number of array of number
        "unit" must be a string
        '''
        if unit == self.reference:
            return value

        for i,known in enumerate(self.accepted):
            if isinstance(known, list):
                for loc_known in known:
                    if unit == loc_known:
                        return self.factor[i]*value + self.offset[i]

            else:
                if unit == known:
                    return self.factor[i]*value + self.offset[i]

        # If unit does not match any accepted unit:
        raise UnknownUnitError('Not recognized unit "'+str(unit)+'"')




#################################################



temperature_units = units('degrees_celsius',
                          [['C', '°', '°C', 'deg C', 'degC', 'deg Celsius', 'deg celsius', 'degree Celsius', 'degree celsius', 'degrees Celsius',
                            'degrees celsius', 'degree_Celsius', 'degree_celsius', 'degrees_Celsius'],
                           ['K', 'Kelvin', 'kelvin']],
                          [1, 1],
                          [0, -273.15])

length_units = units('m',
                     [['meter', 'meters'],
                      ['cm', 'centimeter', 'centimeters'],
                      ['mm', 'millimeter', 'millimeters'],
                      ['km', 'kilometer', 'kilometers'],
                      ['1e6m', '1e6 m', 'Mm', 'megameter', 'megameters'],
                      ['1e6km', '1e6 km', '1e9m', '1e9 m', 'Gm', 'gigameter', 'gigameters'],
                      ['1e9km', '1e9 km', '1e12m', '1e12 m', 'Tm', 'terameter', 'terameters']],
                     [1, 1e-2, 1e-3, 1e3, 1e6, 1e9, 1e12])

area_units = units('m2',
                   [['m^2', 'm**2', 'meter square', 'meter squared', 'square meter', 'squared meter', 'meters square', 'meters squared', 'square meters', 'squared meters'],
                    ['cm^2', 'cm**2', 'centimeter square', 'centimeter squared', 'square centimeter', 'squared centimeter', 'centimeters square', 'centimeters squared', 'square centimeters', 'squared centimeters'],
                    ['mm^2', 'mm**2', 'millimeter square', 'millimeter squared', 'square millimeter', 'squared millimeter', 'millimeters square', 'millimeters squared', 'square millimeters', 'squared millimeters'],
                    ['km^2', 'km**2', 'kilometer square', 'kilometer squared', 'square kilometer', 'squared kilometer', 'kilometers square', 'kilometers squared', 'kilosquare meters', 'squared kilometers', '1e6m2', '1e6 m2', 'M m2', '1e6m^2', '1e6 m^2', 'M m^2', '1e6m**2', '1e6 m**2', 'M m**2'],
                    ['1e9m2', '1e9 m2', 'G m2', '1e9m^2', '1e9 m^2', 'G m^2', '1e9m**2', '1e9 m**2', 'G m**2'],
                    ['1e6km2', '1e6 km2', 'Mkm2', 'M km2', '1e6km^2', '1e6 km^2', 'Mkm^2', 'M km^2', '1e6km**2', '1e6 km**2', 'Mkm**2', 'M km**2', '1e12m2', '1e12 m2', 'T m2', '1e12m^2', '1e12 m^2', 'T m^2', '1e12m**2', '1e12 m**2', 'T m**2'],
                    ['1e9km2', '1e9 km2', 'Gkm2', 'G km2', '1e9km^2', '1e9 km^2', 'Gkm^2', 'G km^2', '1e9km**2', '1e9 km**2', 'Gkm**2', 'G km**2', '1e15m2', '1e15 m2', 'P m2', '1e15m^2', '1e15 m^2', 'P m^2', '1e15m**2', '1e15 m**2', 'P m**2'],
                    ['1e12km2', '1e12 km2', 'Tkm2', 'T km2', '1e12km^2', '1e12 km^2', 'Tkm^2', 'T km^2', '1e12km**2', '1e12 km**2', 'Tkm**2', 'T km**2', '1e18m2', '1e18 m2', '1e18m^2', '1e18 m^2', '1e18m**2', '1e18 m**2']],
                   [1, 1e-4, 1e-6, 1e6, 1e9, 1e12, 1e15, 1e18])

volume_units = units('m3',
                     [['m^3', 'm**3', 'meter cube', 'cubic meter', 'meters cube', 'cubic meters'],
                      ['cm^3', 'cm**3', 'centimeter cube', 'cubic centimeter', 'centimeters cube', 'cubic centimeters'],
                      ['mm^3', 'mm**3', 'millimeter cube', 'cubic millimeter', 'millimeters cube', 'cubic millimeters'],
                      ['1e6m3', '1e6 m3', 'M m3', '1e6m^3', '1e6 m^3', 'M m^3', '1e6m**3', '1e6 m**3', 'M m**3'],
                      ['km^3', 'km**3', 'kilometer cube', 'cubic kilometer', 'kilometers cube', 'cubic kilometers', '1e9m3', '1e9 m3', 'G m3', '1e9m^3', '1e9 m^3', 'G m^3', '1e9m**3', '1e9 m**3', 'G m**3'],
                      ['1e12m3', '1e12 m3', 'T m3', '1e12m^3', '1e12 m^3', 'T m^3', '1e12m**3', '1e12 m**3', 'T m**3'],
                      ['1e6km3', '1e6 km3', 'Mkm3', 'M km3', '1e6km^3', '1e6 km^3', 'Mkm^3', 'M km^3', '1e6km**3', '1e6 km**3', 'Mkm**3', 'M km**3', '1e15m3', '1e15 m3', 'P m3', '1e15m^3', '1e15 m^3', 'P m^3', '1e15m**3', '1e15 m**3', 'P m**3'],
                      ['1e9km3', '1e9 km3', 'Gkm3', 'G km3', '1e9km^3', '1e9 km^3', 'Gkm^3', 'G km^3', '1e9km**3', '1e9 km**3', 'Gkm**3', 'G km**3', '1e18m3', '1e18 m3', '1e18m^3', '1e18 m^3', '1e18m**3', '1e18 m**3']],
                     [1, 1e-6, 1e-9, 1e6, 1e9, 1e12, 1e15, 1e18])

latitude_units = units('degrees_north',
                       [['degree_north', 'degrees_North', 'degree_North', 'degrees north', 'degree north', 'degrees North', 'degree North',
                         'deg N', 'deg North', 'deg north', 'deg', 'degree', 'degrees', '°'],
                        ['rad', 'radian', 'radians']],
                       [1, 3.14159265358979/180])

velocity_units = units('m/s',
                       [['m s^-1', 'm.s^-1', 'm*s^-1', 'm s**-1', 'm.s**-1', 'm*s**-1', 'meter per second', 'meters per second'],
                        ['m/h', 'm h^-1', 'm.h^-1', 'm*h^-1', 'm h**-1', 'm.h**-1', 'm*h**-1', 'm/hr', 'm hr^-1', 'm.hr^-1', 'm*hr^-1', 'm hr**-1', 'm.hr**-1', 'm*hr**-1', 'meter per hour', 'meters per hour'],
                        ['m/d', 'm d^-1', 'm.d^-1', 'm*d^-1', 'm d**-1', 'm.d**-1', 'm*d**-1', 'meter per day', 'meters per day'],
                        ['cm/s', 'cm s^-1', 'cm.s^-1', 'cm*s^-1', 'cm s**-1', 'cm.s**-1', 'cm*s**-1', 'centimeter per second', 'centimeters per second'],
                        ['cm/h', 'cm h^-1', 'cm.h^-1', 'cm*h^-1', 'cm h**-1', 'cm.h**-1', 'cm*h**-1', 'cm/hr', 'cm hr^-1', 'cm.hr^-1', 'cm*hr^-1', 'cm hr**-1', 'cm.hr**-1', 'cm*hr**-1', 'centimeter per hour', 'centimeters per hour'],
                        ['cm/d', 'cm d^-1', 'cm.d^-1', 'cm*d^-1', 'cm d**-1', 'cm.d**-1', 'cm*d**-1', 'centimeter per day', 'centimeters per day'],
                        ['mm/s', 'mm s^-1', 'mm.s^-1', 'mm*s^-1', 'mm s**-1', 'mm.s**-1', 'mm*s**-1', 'millimeter per second', 'millimeters per second'],
                        ['mm/h', 'mm h^-1', 'mm.h^-1', 'mm*h^-1', 'mm h**-1', 'mm.h**-1', 'mm*h**-1', 'mm/hr', 'mm hr^-1', 'mm.hr^-1', 'mm*hr^-1', 'mm hr**-1', 'mm.hr**-1', 'mm*hr**-1', 'millimeter per hour', 'millimeters per hour'],
                        ['mm/d', 'mm d^-1', 'mm.d^-1', 'mm*d^-1', 'mm d**-1', 'mm.d**-1', 'mm*d**-1', 'millimeter per day', 'millimeters per day']],
                       [1, 1/3600, 1/86400, 1e-2, 1e-2/3600, 1e-2/86400, 1e-3, 1e-3/3600, 1e-3/86400])

flux_units = units('m3/s',
                   [['m^3 s^-1', 'm^3/s', 'm^3.s^-1', 'm^3*s^-1', 'm**3/s', 'm**3 s**-1', 'm**3.s**-1', 'm**3*s**-1', 'cubic meter per second', 'cubic meters per second'],
                    ['m3/h', 'm^3/h', 'm^3 h^-1', 'm^3.h^-1', 'm^3*h^-1', 'm**3/h', 'm**3 h**-1', 'm**3.h**-1', 'm**3*h**-1', 'm3/hr', 'm^3/h', 'm^3 hr^-1', 'm^3.hr^-1', 'm^3*hr^-1', 'm**3/hr', 'm**3 hr**-1', 'm**3.hr**-1', 'm**3*hr**-1', 'cubic meter per hour', 'cubic meters per hour'],
                    ['m3/d', 'm^3/d', 'm^3 d^-1', 'm^3.d^-1', 'm^3*d^-1', 'm**3/d', 'm**3 d**-1', 'm**3.d**-1', 'm**3*d**-1', 'cubic meter per day', 'cubic meters per day'],
                    ['Sv', 'sv', 'M m3/s', '1e6m3/s', '1e6 m3/s', '1e6*m3/s', '1e6 m^3 s^-1', 'M m^3 s^-1', 'm^3 s^-1', '1e6m^3/s', '1e6 m^3/s', '1e6*m^3/s', 'M m^3/s', '1e6 m^3.s^-1', '1e6.m^3.s^-1', 'M m^3.s^-1', '1e6 m^3*s^-1', '1e6*m^3*s^-1', 'M m^3*s^-1', '1e6 m**3/s', '1e6m**3/s', '1e6*m**3/s', 'M m**3/s', '1e6 m**3 s**-1', 'M m**3 s**-1', '1e6.m**3.s**-1', 'M m**3.s**-1', '1e6 m**3.s**-1', '1e6 m**3*s**-1', '1e6*m**3*s**-1', 'Sverdrup', 'sverdrup', 'million cubic meter per second', 'million cubic meters per second']],
                   [1, 1/3600, 1/86400, 1e-6])




