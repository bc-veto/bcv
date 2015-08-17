import os
from distutils import log
from distutils.command import build_py
from distutils.command import sdist
from distutils.command import install
from distutils.core import setup

from misc import determine_git_version

def remove_root(path, root):
	if root:
		return os.path.normpath(path).replace(os.path.normpath(root), "")
	return os.path.normpath(path)

class bcv_build_py(build_py.build_py):
  def run(self):
    # create the git_version module
    if determine_git_version.in_git_repository():
      try:
        log.info("generating bcv/git_version.py")
        git_version_fileobj = open("bcv/git_version.py", "w")
        determine_git_version.write_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()
    elif os.path.exists("bcv/git_version.py"):
      # We're probably being built from a release tarball; don't overwrite
      log.info("not in git checkout; using existing bcv/git_version.py")
    else:
      log.info("not in git checkout; writing empty bcv/git_version.py")
      try:
        git_version_fileobj = open("bcv/git_version.py", "w")
        determine_git_version.write_empty_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()

    # resume normal build procedure
    build_py.build_py.run(self)

class bcv_install(install.install):
	def run(self):
		# create the user env scripts
		if self.install_purelib == self.install_platlib:
			bcv_pythonpath = self.install_purelib
		else:
			bcv_pythonpath = self.install_platlib + ":" + self.install_purelib

		bcv_prefix = remove_root(self.prefix, self.root)
		bcv_install_scripts = remove_root(self.install_scripts, self.root)
		bcv_pythonpath = remove_root(bcv_pythonpath, self.root)
		bcv_install_platlib = remove_root(self.install_platlib, self.root)

		if not os.path.isdir("etc"):
			os.mkdir("etc")
		log.info("creating bcv-user-env.sh script")
		env_file = open(os.path.join("etc", "bcv-user-env.sh"), "w")
		print >> env_file, "# Source this file to access BCV"
		print >> env_file, "BCV_PREFIX=" + bcv_prefix
		print >> env_file, "export BCV_PREFIX"
		print >> env_file, "PATH=" + bcv_install_scripts + ":${PATH}"
		print >> env_file, "PYTHONPATH=" + bcv_pythonpath + ":${PYTHONPATH}"
		print >> env_file, "LD_LIBRARY_PATH=" + bcv_install_platlib + ":${LD_LIBRARY_PATH}"
		print >> env_file, "DYLD_LIBRARY_PATH=" + bcv_install_platlib + ":${DYLD_LIBRARY_PATH}"
		print >> env_file, "export PATH PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH"
		env_file.close()

		# now run the installer
		install.install.run(self)


class bcv_sdist(sdist.sdist):
  def run(self):
    # create the git_version module
    if determine_git_version.in_git_repository():
      log.info("generating bcv/git_version.py")
      try:
        git_version_fileobj = open("bcv/git_version.py", "w")
        determine_git_version.write_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()
    else:
      log.info("not in git checkout; writing empty bcv/git_version.py")
      try:
        git_version_fileobj = open("bcv/git_version.py", "w")
        determine_git_version.write_empty_git_version(git_version_fileobj)
      finally:
        git_version_fileobj.close()

    # now run sdist
    sdist.sdist.run(self)

setup(name="BCV",
      version="1.0",
      description="Bilinear Coupling Veto",
      author="Tomoki Isogai, A. Ajith",
      last_modified_by="Bernard Hall, Nairwita Mazumder"
      author_email="isogait@mit.edu, bernard.hall@wsu.edu ",
      packages=["bcv"],
      cmdclass = {
        "build_py": bcv_build_py,
        "install": bcv_install,
        "sdist": bcv_sdist
      },
      scripts=["bin/bcvSetup_unv","bin/bcvOmegaveto","bin/bcvReport","bin/bcvSummaryPage","bin/bcvInsert","bin/bcvWeeklyRun","bin/getKW","compile/omegaveto"],
      data_files = [("etc",["etc/bcv-user-env.sh"]),
                    ("config",["config/example.ini","config/configuration.txt","config/weekly.ini"])]
     )
