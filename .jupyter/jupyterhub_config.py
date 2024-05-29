
c = get_config()

def email_lookup_hook(spawner, auth_state):
    # TODO: once the firstuseauthenticator is replaced by ldapauthenticator,
    # the whole function should simplify to something like this:
    # spawner.environment['NXF_MAIL'] = auth_state['mail']

    import ldap
    LDAP_URL = "ldaps://ldap-01.zdv.uni-mainz.de"
    LDAP_DN = "dc=uni-mainz,dc=de"
    LDAP_SCOPE = ldap.SCOPE_SUBTREE
    LDAP_FILTER = f'(uid={spawner.user.name})'
    
    l = ldap.initialize(LDAP_URL)
    l.set_option(ldap.OPT_X_TLS_REQUIRE_CERT, 0)
    l.set_option(ldap.OPT_X_TLS_NEWCTX, 0)
    
    spawner.environment['NXF_MAIL'] = l.search_s(LDAP_DN, LDAP_SCOPE, LDAP_FILTER)[0][1]['mail'][0].decode('utf-8')

# authenticator config
c.JupyterHub.authenticator_class = 'firstuseauthenticator.FirstUseAuthenticator'
c.Authenticator.allow_all = True

# network config
c.JupyterHub.bind_url = 'http://:42420'
c.ConfigurableHTTPProxy.api_url = 'http://localhost:42421'
c.JupyterHub.hub_port = 42424
c.JupyterHub.hub_ip = 'hpc1'

# spawner config
c.JupyterHub.spawner_class = 'nextflow'
c.Spawner.auth_state_hook = email_lookup_hook
c.Spawner.args = ['--debug']
c.Spawner.debug = True
c.Spawner.env_keep = ['PATH', 'PYTHONPATH', 'CONDA_ROOT', 'CONDA_DEFAULT_ENV', 'VIRTUAL_ENV', 'LANG', 'LC_ALL', 'JUPYTERHUB_SINGLEUSER_APP', 'PIXI_PROJECT_MANIFEST', 'PIXI_ENVIRONMENT_NAME', 'PIXI_PROJECT_NAME', 'PIXI_PROJECT_ROOT']
c.Spawner.http_timeout = 600

