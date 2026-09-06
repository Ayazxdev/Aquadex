import React, { useState, useEffect } from 'react';
import { Beaker, Home, FlaskConical, User, LogOut } from 'lucide-react';
import Login from '../components/Login'; // Adjust path as needed
import './header.css';

export default function Header({ currentPage, onNavigate }) {
  const [isScrolled, setIsScrolled] = useState(false);
  const [showLogin, setShowLogin] = useState(false);
  const [user, setUser] = useState(null);

  useEffect(() => {
    const handleScroll = () => {
      setIsScrolled(window.scrollY > 20);
    };

    window.addEventListener('scroll', handleScroll);
    
    // Check for existing user session
    const savedUser = localStorage.getItem('currentUser');
    if (savedUser) {
      setUser(JSON.parse(savedUser));
    }

    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  const handleLogin = (userData) => {
    setUser(userData);
    setShowLogin(false);
  };

  const handleLogout = () => {
    setUser(null);
    localStorage.removeItem('currentUser');
    localStorage.removeItem('rememberedUser');
  };

  return (
    <>
      <header className={`header ${isScrolled ? 'scrolled' : ''}`}>
        <div className="header-container">
          <div className="header-left">
            <div className="header-logo" onClick={() => onNavigate('home')}>
              <Beaker className="logo-icon" />
              <span className="logo-text">AQUADEX</span>
            </div>
            
            <nav className="header-nav">
              <button 
                className={`nav-link ${currentPage === 'home' ? 'active' : ''}`}
                onClick={() => onNavigate('home')}
              >
                Home
              </button>
              <button 
                className={`nav-link ${currentPage === 'run' ? 'active' : ''}`}
                onClick={() => onNavigate('run')}
              >
                Run
              </button>
              <button 
                className={`nav-link ${currentPage === 'results' ? 'active' : ''}`}
                onClick={() => onNavigate('results')}
              >
                Results
              </button>
            </nav>
          </div>

          <div className="header-right">
            <button 
              className="btn-header-primary"
              onClick={() => onNavigate('run')}
            >
              Start Analysis
            </button>
            
            {user ? (
              <div className="user-menu">
                <div className="user-info">
                  <User size={16} />
                  <span className="user-name">{user.name}</span>
                </div>
                <button 
                  className="btn-logout"
                  onClick={handleLogout}
                  title="Logout"
                >
                  <LogOut size={16} />
                </button>
              </div>
            ) : (
              <button 
                className="btn-login"
                onClick={() => setShowLogin(true)}
              >
                Login
              </button>
            )}
          </div>
        </div>
      </header>

      {/* Login Modal */}
      {showLogin && (
        <Login 
          onClose={() => setShowLogin(false)} 
          onLogin={handleLogin}
        />
      )}
    </>
  );
}