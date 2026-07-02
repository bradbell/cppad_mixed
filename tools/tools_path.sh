# ---------------------------------------------------------------------------
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
# SPDX-FileContributor: 2026 Bradley M. Bell
# ---------------------------------------------------------------------------
# source tools/tools_path.sh
# Adds the local directory tools, below the current directory, to the PATH
# ---------------------------------------------------------------------------
PATH=$(echo $PATH | sed -e 's|^|tools:|' -e 's|:tools:|:|')
