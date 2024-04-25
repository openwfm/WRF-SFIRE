
The WRF-Fire code in [WRF](https://wrf-model.org/WRF) release has an old version of the code here, from 2010, with further modifications at NCAR and developments since WRF 4.0. We cannot support that version and it is not under our control. Many different features were added since then. Please get the current version from https://github.com/openwfm/WRF-SFIRE as described here, which includes current SFIRE and WRF. This is the only version we can support. We can even consider adding features on request.

SFIRE is maintained in a git repository. A git repository is a software version control system similar to CVS or SVN. There is a repository containing WRF-SFIRE and another repository containing WPS, and it is regularly synchronized with major WRF releases. They both need to be compiled with the same compiler and be present at the same directory level. In order to obtain the software from our repositories, you must first download and install the git software. The default view of the repository is the head of the master branch, which should contain a current stable and working code at any time.

Developers can get write access to the git repository on request, and must push changes to it.

## Download

We strongly recommend using download by git (as opposed to a tar or zip file, also available from the repository), because git allows easy updates and identification of the state of all files in case of problems. The following public, read-only mirrors are available:

- [https://github.com/openwfm/WRF-SFIRE](http://github.com/openwfm/WRF-SFIRE): `git clone git://github.com/openwfm/WRF-SFIRE.git`
- [http://repo.or.cz/WRF-SFIRE.git](http://repo.or.cz/WRF-SFIRE.git): `git clone git://repo.or.cz/WRF-SFIRE.git`

The home page of the mirror allows you to browse the source code as well as examine changes and view the different development branches.

It is also possible to download the latest source code as a tarball or zip file from each mirror; however, this method is discouraged as it will be difficult to update to the latest version and we cannot support code that is not under git control.

## Update

You can update your files at any time to the current version without downloading everything again. If you have downloaded your copy by git, do in your WRF-SFIRE directory:

```bash
git checkout master; git pull
```

If you have changed anything, the update will fail, and you have to commit your changes first by `git commit -a`. The `git pull` will then merge your changes with the new files from the repository. Please be sure you work with the latest files before contacting us for support. However, we may not be able to support code that you have changed.

## Verify

To identify the version of files you have, do in your WRF-SFIRE directory:

```bash
git log
```

To verify your files, use:

```bash
git diff
```

If your files have not changed from the version identified by `git log`, `git diff` will give no output.
