# WRF-Fire, SFIRE, and WRF-SFIRE Lineage

This note summarizes the relationship between WRF-Fire, SFIRE, and WRF-SFIRE as it appears from the OpenWFM wiki, the original `openwfm/wrf-fire` repository, and the branch/tag history preserved in this repository.

## Naming

The names are historically intertwined, but they are useful to distinguish:

- **SFIRE** is the fire spread model/code lineage.
- **WRF-SFIRE** is the OpenWFM-maintained coupled WRF plus SFIRE system.
- **WRF-Fire** usually refers to the fire module distributed inside official WRF releases.

The OpenWFM wiki makes this split explicitly: the actively developed fire code is called SFIRE, the coupled system is WRF-SFIRE, and the name WRF-Fire is now generally reserved for the WRF-release version.

## Historical Arc

The earliest fire history in this repository points to the WRF plus fire code entering the tree in 2007:

```text
2007-06-06  b7ed4f2e4  initial version of wrf+fire code given by Ned Patton at NCAR.
```

By late 2007, the history already contains `sfire_driver`, `sfire_model`, WRF coupling, `em_fire`, tile work, and Matlab-oriented standalone tests. This matches the OpenWFM history: WRF-Fire coupling work from CAWFE/NCAR was developed into the SFIRE level-set implementation and its WRF integration.

A key handoff to official WRF is preserved as:

```text
2010-01-15  WRF-Fire-included-in-3.3
bbf025b7b   Updates to the WRF Fire code, provided by Jonathan Beezely,
            Jan Mandel, and Janice Coen. This version includes modifications
            by John Michalakes to permit subgrid history, conforming to WRF
            best practices.
```

That commit changes `Registry/registry.fire`, `phys/module_fr_sfire_*`, and `test/em_fire`, so the code handed into WRF 3.3 was still internally using the SFIRE names.

After that, two lines continued:

- official WRF kept a WRF-Fire module in WRF releases;
- OpenWFM kept developing the SFIRE/WRF-SFIRE line separately.

There is also a revealing rename point:

```text
2013-01-20  renamed-sfire-to-fire
ff7ec687b   Rename all of the fire files, and the associated subroutines.
```

That commit renames files such as:

```text
phys/module_fr_sfire_driver.F -> phys/module_fr_fire_driver.F
phys/module_fr_sfire_phys.F   -> phys/module_fr_fire_phys.F
```

In other words, one WRF-release-oriented line moved the internal names back toward `fire`, while the OpenWFM development line retained SFIRE as the name of the actively developed fire model.

## Repository Roles

The original historical repository still exists:

- `openwfm/wrf-fire`: <https://github.com/openwfm/wrf-fire>

Its README describes it as the original WRF-Fire repository and says development moved to:

- `openwfm/WRF-SFIRE`: <https://github.com/openwfm/WRF-SFIRE>

The original repo is useful for archaeology. It preserves the older layout, including `wrfv2_fire`, `WRFV3`, `standalone`, `WPS`, and documentation. The current WRF-SFIRE repository is the maintained WRF fork that carries the evolving SFIRE code.

## Branch Head Map

The current local repository contains several imported/tracking histories. At the time this was checked, the useful heads were:

```text
origin/WRF-track/master
  a8eb84685  2023-12-21
  Merge remote-tracking branch 'origin/release-v4.5.2'

origin/wrf-fire-track/master
  a14489984  2019-12-02
  adding configure.wrf.intel.debug

origin/wrf-fire-track/submitted-to-3.3
  3cdb542af  2019-01-02
  fixed error in heat flux instertion code that was there from the beginning

origin/master
  9801186f3  2026-05-01
  Fast forward to develop-3 standalone. Close #3.

origin/release_S02
  7b5376646  2025-10-27
  remove no_fuel_cat and no_fuel_cat2
```

The main historical tags/milestones are:

```text
WRF-Fire-included-in-3.3  bbf025b7b  2010-01-15
renamed-sfire-to-fire     ff7ec687b  2013-01-20
v4.0-fire                 e26f21fb3  2018-03-09
submitted-to-3.3          3cdb542af  2019-01-02
wrf-sfire-v4.3.1          7a0e87f7e  2021-11-06
merge-v4.4                4c9bccda8  2022-05-11
W4.4-S0.1                 c95958009  2024-01-12
```

## Rough Divergence

These counts are branch-head comparisons, not semantic summaries:

```text
origin/WRF-track/master...origin/master
  96 commits only on WRF-track/master
  1877 commits only on WRF-SFIRE master

origin/wrf-fire-track/master...origin/master
  2714 commits only on wrf-fire-track/master
  8333 commits only on WRF-SFIRE master

origin/wrf-fire-track/submitted-to-3.3...origin/master
  1497 commits only on wrf-fire submitted-to-3.3 side
  8333 commits only on WRF-SFIRE master

origin/master...origin/release_S02
  14 commits only on master
  89 commits only on release_S02
```

The high-level conclusion is that current WRF-SFIRE is not simply "WRF plus original wrf-fire." It is a long-running WRF fork that periodically imports WRF baselines while also carrying a separately evolved SFIRE line. The original `openwfm/wrf-fire` history survives mostly as historical/tracking/filter branches, useful for archaeology and patch comparison rather than as the active trunk.

## Sources

- OpenWFM WRF-SFIRE wiki: <https://wiki.openwfm.org/wiki/WRF-SFIRE>
- OpenWFM WRF-Fire wiki: <https://wiki.openwfm.org/wiki/WRF-Fire>
- How to get WRF-SFIRE: <https://wiki.openwfm.org/wiki/How_to_get_WRF-SFIRE>
- Original repository: <https://github.com/openwfm/wrf-fire>
- Current repository: <https://github.com/openwfm/WRF-SFIRE>
- Official WRF-Fire user guide: <https://www2.mmm.ucar.edu/wrf/site/documentation/users_guide/fire.html>

## Provenance

Generated automatically from the branch, tag, and commit history in this repository, with supporting context from the OpenWFM wiki and linked upstream repositories.

Signed-off-by: Codex <codex@openai.com>
