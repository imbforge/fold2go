c = get_config()

# authenticator config
c.JupyterHub.authenticator_class = 'firstuseauthenticator.FirstUseAuthenticator'

# network config
c.JupyterHub.bind_url = 'http://:42420'
c.ConfigurableHTTPProxy.api_url = 'http://localhost:42421'
c.JupyterHub.hub_port = 42424
c.JupyterHub.hub_ip = 'hpcgpu'

# spawner config
c.JupyterHub.spawner_class = 'nextflow'
c.Spawner.args = ['--debug']
c.Spawner.debug = True
c.Spawner.env_keep = ['PATH', 'PYTHONPATH', 'CONDA_ROOT', 'CONDA_DEFAULT_ENV', 'VIRTUAL_ENV', 'LANG', 'LC_ALL', 'JUPYTERHUB_SINGLEUSER_APP', 'PIXI_PROJECT_MANIFEST', 'PIXI_ENVIRONMENT_NAME', 'PIXI_PROJECT_NAME', 'PIXI_PROJECT_ROOT']
c.Spawner.http_timeout = 60