# AIRES `geomagnetic_()` interface evidence

AIRES is outside this repository. All paths below are relative to `$REPO/..`, where `$REPO` is `MaGICS.1-5-0`.

## Implementation location

The source implementation of `geomagnetic_()` was not present in the readable AIRES tree. Evidence found:

| Evidence | Path | Lines |
| --- | --- | --- |
| MaGICS declares and calls `geomagnetic_(&theta,&phi,&altitude,&fdate,&b,&inc,&dec)` | `MaGICS.1-5-0/src/MaGICS/geometry.c` | 67-72 |
| AIRES static library exports `geomagnetic_` | `aires/19-04-08/lib/libAires.a` | `nm`: symbol in `AiresCalls.o` |
| AIRES shared library exports `geomagnetic_` | `aires/19-04-08/lib/libAires.so` | `nm -D`: symbol at `00000000000800d0` |
| Disassembly of `geomagnetic_()` calls `cvtfltdate_` and `igrf13plsyn_` | `aires/19-04-08/lib/libAires.so` | `objdump` around symbol `00000000000800d0` |

The binary evidence shows that the routine converts the floating date and then calls the IGRF synthesis routine. The manual says AIRES evaluates the magnetic field with IGRF when geographic place and date are specified.

## Inferred C/Fortran signature

From the MaGICS call site, the interface used here is:

```c
void geomagnetic_(double *theta,
                  double *phi,
                  double *altitude,
                  double *fdate,
                  double *b,
                  double *inclination,
                  double *declination);
```

Inputs:

| Argument | Meaning used by MaGICS | Unit |
| --- | --- | --- |
| `theta` | geographic latitude | degrees |
| `phi` | geographic longitude in AIRES convention | degrees |
| `altitude` | altitude above sea level | meters |
| `fdate` | event date as AIRES floating date | year or converted `year month day` |

Outputs:

| Argument | Meaning | Unit |
| --- | --- | --- |
| `b` | geomagnetic field strength `F = |B|` | nanotesla |
| `inclination` | magnetic inclination `I` | degrees |
| `declination` | magnetic declination `D` | degrees |

The first diagnostic run shows that `b` is returned numerically in nT: AIRES reports `22.665 uT` near the site, while the MaGICS diagnostic CSV reports `B = 22661.8` at the same point; `22661.8 nT = 22.6618 uT`. This matches the public AIRES input-data convention for magnetic field strength.

Important consequence: MaGICS' `chi()` uses this nT-valued `Bperp` directly against `bcrit`, whose numeric value is computed from SI constants at `MaGICS.1-5-0/src/MaGICS/init.c:169-170`. The current `chi()` remains numerically correct because `particle.energy` is numerically in GeV and the nT magnetic field contributes the reciprocal factor `1e-9`.

## Units of `B`

Manual evidence:

| Evidence | Path | Lines |
| --- | --- | --- |
| Magnetic field is described by strength `F = |B|`, inclination `I`, and declination `D` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 1293-1297 |
| Accepted magnetic units include tesla, gauss, and gamma/nT | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 2989-2993 |
| Manual `GeomagneticField 32 uT -60 deg 2 deg` example maps parameters to `F`, `I`, `D` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4114-4127 |
| Public AIRES `realdata[22]` is field strength in nT | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 9461-9465 |
| MaGICS diagnostic CSV gives `B = 22661.8` where AIRES reports `22.665 uT` | first diagnostic run | runtime evidence |

Conclusion for MaGICS: `geomagnetic_()` output `b` is numerically in nT. Any future unit cleanup must convert both the photon energy and magnetic field conventions together, because changing only one side would change `chi` by `1e9`.

## Units of `altitude`

Manual evidence:

| Evidence | Path | Lines |
| --- | --- | --- |
| Site altitude is altitude above sea level | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4021-4024 |
| Predefined site table labels altitude as `m.a.s.l` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4080-4097 |
| Public AIRES input data array lists injection, ground and observing altitudes in meters | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 9451-9458 |

MaGICS passes `altitude` in meters: positions and `earth_radius` are in meters (`MaGICS.1-5-0/src/MaGICS/constants.h:3-10`), `chi()` computes altitude from `earth_radius` at `MaGICS.1-5-0/src/MaGICS/geometry.c:166`, and then passes it to `cartesian_magnetic_field()` at `MaGICS.1-5-0/src/MaGICS/geometry.c:168-171`.

## Conventions for `theta` and `phi`

MaGICS local spherical coordinates are documented at `MaGICS.1-5-0/src/MaGICS/geometry.c:21-33`:

| Symbol | MaGICS meaning |
| --- | --- |
| `theta` | latitude, positive in the northern hemisphere |
| `phi` | longitude, positive west |

The same comment warns that MaGICS longitude is opposite to the AIRES convention. Accordingly, `cartesian_magnetic_field()` does:

```c
theta=latitude;
phi=-longitude;
geomagnetic_(&theta,&phi,&altitude,&fdate,&b,&inc,&dec);
```

at `MaGICS.1-5-0/src/MaGICS/geometry.c:67-72`.

AIRES manual evidence for longitude sign appears in the site table, where longitudes are reported as east or west (`aires/19-04-08/doc/UsersManual190408.pdf` text extraction lines 4064-4078), and the example `AddSite cld -31.5 deg -64.2 deg 387 m` lists latitude, longitude and altitude at lines 4021-4024. This supports AIRES' signed geographic longitude convention, with west negative in that example.

## Conventions for inclination and declination

Manual evidence:

| Quantity | Convention | Path | Lines |
| --- | --- | --- | --- |
| Inclination `I` | angle between local horizontal plane and field vector; positive when `B` points downwards | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 1293-1297 |
| Declination `D` | angle between horizontal field component and geographical north; positive toward east | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 1293-1297 |
| AIRES x-axis | local magnetic north, direction of horizontal component of geomagnetic field | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 1308-1309 |
| Manual geomagnetic directive parameters | `F`, `I`, `D` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 7384-7388 |

MaGICS uses `inc` and `dec` in degrees, converting with `pi/180` at `MaGICS.1-5-0/src/MaGICS/geometry.c:73-76`.

## Format and meaning of `fdate`

Manual evidence:

| Evidence | Path | Lines |
| --- | --- | --- |
| `Date` syntax accepts either floating year or `year month day` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4025-4028 |
| Floating date uses year as unit; second format is `year month day` | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4098-4099 |
| Date is used to evaluate the geomagnetic field with IGRF | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 7141-7151 |
| If date is missing, AIRES uses system time at simulation start | `aires/19-04-08/doc/UsersManual190408.pdf` text extraction | 4105-4108 |
| Example: `SetGlobal EventDate 2018 12 31`; `Date {EventDate}` | `aires/19-04-08/demos/01RunningAires/04MultiPrimary/MyEvent001` | 23-31 |
| Example: `Date 2019 04 15` | `aires/19-04-08/demos/01RunningAires/04MultiPrimary/samplevt002.inp` | 28-30 |

MaGICS reads `fdate = realdata[26]` at `MaGICS.1-5-0/src/MaGICS/init.c:150`. The disassembly of `geomagnetic_()` shows it calls `cvtfltdate_`, consistent with passing the AIRES floating date to the low-level routine.
