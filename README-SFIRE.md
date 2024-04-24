# WRF-SFIRE
This repository is https://github.com/openwfm/WRF-SFIRE master branch. 

A mirror with a graphical log is at https://repo.or.cz/git-browser/by-commit.html?r=WRF-SFIRE.git


# Table of Contents
- [Current Version](#current-version)
- [Introducing Semantic Versioning](#introducing-semantic-versioning)
- [Branch Management](#branch-management)
- [Regression Testing Introduction](#regression-testing-introduction)
- [Code Quality](#code-quality)
- [Version 1.0: Feature Highlights](#version-10-feature-highlights)
- [Issues](#issues)
- [Branches](#branches)
- [How to run](#how-to-run)
- [Notes](#notes)
- [Notes for developers](#notes-dev)

  
# Current version <a name="current-version"></a>

The current released version of **WRF-SFIRE** is the version `W4.4-S0.1`.
*Release information about this version is not available yet*.
The used version of WRF is 4.4 (see [release information](https://github.com/wrf-model/WRF/releases/tag/v4.4)).
The version of SFIRE is 0.1.

# Introducing Semantic Versioning: Transitioning from Legacy Code to Version 0.1 <a name="introducing-semantic-versioning"></a>

## Background
We are excited to share a significant development in our project's versioning system. As we release the official 0.1 version of our code, we want to introduce you to an enhanced versioning approach that incorporates the relationship between SFIRE code and WRF code. Moving forward, our versioning system will consist of two components: the Weather Research and Forecasting (WRF) code version and the Wildfire Spread Simulation (SFIRE) code version.

## Understanding the Versioning Scheme
The semantic versioning scheme (https://semver.org) we have adopted is designed to provide clarity and context about the interconnectedness of our code with the WRF and SFIRE codes. Each version will be represented in the following format: `Wx.y-Si.j`, where:

- `Wx.y` denotes the version of the WRF code on which our codebase relies. For example, `W4.4` represents version 4.4 of the WRF code. This component indicates the underlying framework that drives our code's functionality.

- `Si.j` represents the version of our code. For instance, `S0.1` signifies the first official release (version 0.1) of our code, SFIRE. This component focuses on the development and evolution of our code specifically.

## Benefits of the Enhanced Versioning System 
This enhanced versioning system offers several advantages to our community:

1. **Clear relationship between codes**: By incorporating the WRF version alongside our code's version, we provide a transparent representation of the underlying codebase that supports our functionality. This clarity ensures that users can easily identify the connection between different versions.

2. **Compatibility and integration**: The WRF-SFIRE integration necessitates version compatibility between the two codes. By explicitly denoting the WRF version, we ensure that users can determine which versions of WRF are compatible with specific releases of our code, thereby facilitating seamless integration.

3. **Streamlined updates**: With the enhanced versioning system, it becomes simpler to track and understand the evolution of both the WRF and SFIRE codes. This knowledge empowers users to make informed decisions about updating either component, based on their specific requirements and compatibility constraints.

4. **Precise communication**: When discussing our code or reporting issues, the enhanced versioning system enables clear and concise communication by explicitly referring to the WRF and SFIRE versions being utilized. This precision helps in resolving potential issues more efficiently.

# Branch Management <a name="branch-management"></a>

## Introduction 
We have implemented a new branch management strategy to improve code organization and streamline the development process. This strategy involves the use of specific branches with defined purposes and guidelines. The key branches in our repository are:

- **Master Branch:** The `master` branch contains officially released content and serves as a stable version of the codebase. Commits on this branch are tagged with specific release tags following the standard format: `Wx.y-Si.j`.

- **Release Branch:** The `release` branch is used to merge developments of SFIRE from `develop-x` branches. It acts as an intermediary branch before merging into the `master` branch. Merging the `release` branch into `master` leads to an increase in the `Si.j` part of the version number.

- **Merge-WRF Branch:** The `merge-WRF` branch is responsible for merging upgrades made on the WRF side to the WRF-SFIRE codebase. This branch facilitates the integration of improvements from the WRF repository into our codebase. Merging the `merge-WRF` branch into `master` will increase the `Wx.y` part of the version number.

- **Develop-x Branches:** The `develop-x` branches are issue-based branches, where `x` refers to an open issue on GitHub. Each `develop-x` branch focuses on addressing a specific issue. Once a `develop-x` branch has been merged into the `release` branch and subsequently into the `master` branch, the related issue is closed automatically by the corresponding commit.

## Branch Workflow - Upgrading SFIRE
The branch workflow for upgrading SFIRE can be summarized as follows:

1. **Feature Development:** Developers create issue-based branches (`develop-x`) to work on specific issues identified on GitHub related to SFIRE enhancements.

2. **Merge to Release Branch:** Once the development work on an issue is completed, the corresponding `develop-x` branch is merged into the `release` branch.

3. **Release to Master Branch:** Periodically, the changes from the `release` branch are merged into the `master` branch to create stable releases. This merge increases the `Si.j` version number, where `i` represents the SFIRE version.

By following this workflow, we ensure that SFIRE enhancements undergo thorough development, testing, and integration before being released in a stable version.

## Branch Workflow - Upgrading WRF Version
The branch workflow for upgrading the WRF version can be summarized as follows:

1. **Merge WRF Upgrades:** Updates from the WRF repository are merged into the `merge-WRF` branch to incorporate improvements from the WRF side into our codebase.

2. **Merge to Master Branch:** The `merge-WRF` branch is periodically merged into the `master` branch to integrate the WRF upgrades. This merge increases the `Wx.y` version number, where `x` represents the major WRF version and `y` represents the minor WRF version.

By following this workflow, we ensure that the WRF-SFIRE codebase remains up to date with the latest enhancements and features from the WRF project.

## Regression Testing

Before any merge into the `master` branch, we apply a rigorous regression testing protocol, described in the [next section](#regression-testing-introduction).

By following this branch management workflow and incorporating regression testing, we aim to ensure a controlled and systematic approach to releasing stable versions while maintaining code integrity.

# Regression Testing Introduction <a name="regression-testing-introduction"></a>

## Introduction
In our ongoing commitment to quality assurance and delivering a reliable codebase, we are excited to introduce regression testing as an integral part of our development process. Regression testing is a systematic approach to retesting previously implemented functionalities to ensure that recent changes or additions to the codebase have not introduced unintended side effects or regressions.

## Benefits of Regression Testing
By incorporating regression testing, we aim to achieve the following benefits:

1. **Ensure Code Stability:** Regression testing allows us to verify that existing functionalities continue to work as intended after implementing new features or making code modifications. This helps us identify and fix any unforeseen issues that may have arisen due to recent changes.

2. **Prevent Unintended Side Effects:** Through a comprehensive suite of regression tests, we can proactively detect any unintended side effects or regressions that could potentially impact the performance, reliability, or functionality of our codebase.

3. **Safeguard Against Future Changes:** By maintaining a robust regression test suite, we establish a safety net that protects against future code modifications or updates. The tests serve as a baseline to ensure that future changes do not inadvertently break existing functionalities.

## Implementation Approach
To implement regression testing, our development team will:

- Create a suite of well-defined and comprehensive tests that cover critical functionalities, edge cases, and frequently used features.
- Automate the execution of these tests to streamline the testing process, improve efficiency, and allow for continuous integration.
- Regularly update and expand the regression test suite to encompass new features, bug fixes, and code enhancements introduced in each release.

# Code Quality <a name="code-quality"></a>

## Introduction
In our ongoing efforts to maintain a high-quality codebase, we are introducing coding and linting rules that are aligned with the WRF (Weather Research and Forecasting) coding rules. These rules ensure consistent and reliable code across our project, following industry best practices and the established standards of the WRF community.

## Coding and Linting Rules
We have established a set of coding and linting rules that are specifically designed to comply with the [WRF coding rules](https://www2.mmm.ucar.edu/wrf/WG2/WRF_conventions.html) and inpired by [FORTRAN90/95 coding conventions](https://alm.engr.colostate.edu/cb/wiki/16983). These rules cover various aspects of coding style, naming conventions, documentation practices, and more. By adhering to these rules, we aim to improve code consistency, readability, and collaboration within the project, while also maintaining compatibility with the broader WRF ecosystem.

To enforce these rules, we have integrated automated code analysis tools into our development process. These tools will check the code against the defined rules and provide immediate feedback to developers, allowing them to address any violations and ensure that the codebase aligns with the WRF coding standards.

## Regular Code Quality Evaluation
In addition to the introduction of coding and linting rules, we are implementing a regular code quality evaluation process. This process involves periodic assessments of the codebase to identify areas for improvement, detect potential issues, and ensure compliance with the WRF coding guidelines. To enforce coding standards and best practices, we utilize the [Flint Python library](https://gitlab.com/cerfacs/flint). Flint allows us to perform automated code analysis and evaluate code quality quantitatively. In addition to the default rules provided by Flint, we have implemented custom rules specific to our project. These rules align with the WRF coding guidelines and reflect our coding conventions.


Through this evaluation process, we will proactively address any code quality concerns, optimize performance, enhance maintainability, and strengthen the overall stability of the codebase. By aligning with the WRF coding rules, we aim to uphold the highest standards of code quality and facilitate seamless integration with the wider WRF community.

## Community Engagement
We highly value the involvement and expertise of our community in maintaining code quality. We encourage all community members to embrace the coding and linting rules based on the WRF coding rules, provide feedback, and contribute to the continuous improvement of our codebase.

Your contributions, such as code reviews, suggestions, and bug reports, are invaluable in helping us ensure the highest level of code quality and deliver an exceptional experience to our users, while maintaining compatibility and coherence with the WRF ecosystem.

# SFIRE version 0.1: Feature Highlights <a name="version-10-feature-highlights"></a>

## Level-Set Fire Spread
- Level set function is propagated with RK2 (Heun's) and ENO1 schemes
- Wind and slope are projected on the normal direction with a 2D method.

## Rate of Spread Parameterizations
- Legacy CAWFE version of [Rothermel rate of spread (ros) model](https://doi.org/10.2737/RMRS-GTR-371).
- Rothermel model adapted from [BEHAVE](https://www.frames.gov/behaveplus/home) model.
- McArthur model.
<!-- - Adapted Balbi model from [Chatelon et al. (2022)](https://doi.org/10.1071/WF21082). -->
<!-- - [Scott ros model](https://doi.org/10.2737/RMRS-RP-29) for crown fires. -->

## Surface Fuel Models
- [13 Anderson](https://doi.org/10.2737/INT-GTR-122) fuel models.
<!-- - [40 Scott and Burgan](https://doi.org/10.2737/RMRS-GTR-153) fuel models. -->

## Ignition Pattern
- Point ignition.
- Line ignition.
- Arrival time ignition.

## Heat Flux Parameterizations
- Surface heat flux with predetermined exponential decay.
<!-- - Canopy heat flux with predermined exponential decay independant from the surface. -->

## Fuel Moisture Parameterizations
- Time lag model with 1h, 10h and 100h fuels.

## Surface Initialization
- Surface fuel model map.
- Canopy fuel data map.

## Other Parameterizations
- [Massman](https://doi.org/10.1139/cjfr-2016-0354) canopy wind parameterization.

## Miscellaneous
- Write fuels data into fuels.m file (development feature)

# Issues <a name="issues"></a>
## Linking issues
Link to issues always using:
organization/repository#number
For instance: `wirc-sjsu/WRF-SFIRE#4` or `openwfm/WRF-SFIRE#86`
This will prevent linking to different issues in both organizations with the same number. Links created for commits dated before the issue seems to not cause a problem.

# Branches <a name="branches"></a>
* master, develop - the merged code.
* WRF-track/master, WRF-track/develop, etc. - tracking the official WRF repository including release tags
* filtered/* - rebased commits from the original wrf-fire repository, directory wrfv2_fire mapped to root:
_  filtered/balbi
-  filtered/devel
-  filtered/dvm_branch
-  filtered/master
-  filtered/submitted-to-3.3
-  filtered/wrf-fire-branch-for-fmc-merge
-  fuel-moisture-model
* wrf-fire-track/* - the original wrf-fire
-  wrf-fire-track/master
-  wrf-fire-track/dvm_branch


# How to run <a name="how-to-run"></a>
* In namelist input, ifire=1 is WRF-SFIRE, ifire=2  the fire code in WRF we put there in 2012 (see branch filtered/submitted-to-3.3) with changes at NCAR. Most namelist flags are the same but each version has some of its own. 
* Test problems available in test/em_fire:
- rain ifire=2
- hill ifire=1
_ two_fires ifire=2
* The problems with ifire=2 should give the same result in branch master (the merged code) and in branch WRF/develop (the WRF distribution)
* It seems to run serial and MPI and so far passed some limited testing on cheyenne.

# Notes <a name="notes"></a>
Standalone is included but have not started on it yet.
Branches balbi and dvm_branch were carried over but they are not merged into master or tested yet because they were not merged into master in wrf-fire. Also I do not know how to test dvm_branch. 
I can't do real problems yet, that requires pushing data through current WPS, which Adam and Angel are working on.

# Notes for developers <a name="notes-dev"></a>
## How to upgrade WRF version 
* git checkout WRF/master
* git pull git@github.com:wrf-model/WRF.git master
* git checkout develop
* git merge WRF/master
* resolve conflicts
* copy README from wrf-model to README-WRF
* Make sure that README.md softlink still points to README-SFIRE.md
* run regression tests https://github.com/openwfm/regression
* git checkout master
* git merge develop
